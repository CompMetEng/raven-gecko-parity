function results = curate_tables_edge_toy(ctx)
% MATLAB side of the curateModelFromTables edge-case scenario.
%
% The model has one reaction and two metabolites with a formula and a charge. Each case curates it from a
% set of tables and records the metabolites (formula and charge), the reactions (bounds and stoichiometry)
% and the genes by name, so that identifiers, which each side numbers its own way, are not compared.

inputs = ctx.inputs;
cases = as_list(inputs.cases);
results = struct();
for k = 1:numel(cases)
    results.(char(cases{k}.name)) = curate(inputs, cases{k});
end
end


function results = curate(inputs, c)
model = struct();
model.id = 'toy';
model.rxns = {'R0'};
model.rxnNames = {'existing'};
model.mets = {'x_c'; 'y_c'};
model.metNames = {'x'; 'y'};
model.metFormulas = {'CH4'; 'C2H6'};
model.metCharges = [1; 2];
model.S = sparse([-1; 1]);
model.lb = 0;
model.ub = 1000;
model.rev = 0;
model.c = 0;
model.b = zeros(2, 1);
model.comps = {'c'};
model.compNames = {'c'};
model.metComps = [1; 1];
model.genes = cell(0, 1);
model.grRules = {''};
model.rxnGeneMat = sparse(1, 0);

folder = char(inputs.tables_dir);
args = {};
names = {'genes', 'genesInfo'; 'rxns', 'rxnsInfo'; 'coeffs', 'rxnsCoeffs'};
for q = 1:size(names, 1)
    name = char(c.(names{q, 1}));
    if ~strcmp(name, 'none')
        args = [args, {names{q, 2}, fullfile(folder, name)}]; %#ok<AGROW>
    end
end
metsInfo = char(c.mets);
if ~strcmp(metsInfo, 'none')
    metsInfo = fullfile(folder, metsInfo);
end
args = [args, {'metPrefix', char(inputs.met_prefix), 'rxnPrefix', char(inputs.rxn_prefix)}];
if isfield(c, 'overwrite_with_empty')
    args = [args, {'overwriteWithEmpty', logical(c.overwrite_with_empty)}];
end

try
    evalc('newModel = curateModelFromTables(model, metsInfo, args{:});');
catch
    results = struct('error', true);
    return
end

labels = cell(numel(newModel.mets), 1);
for i = 1:numel(labels)
    labels{i} = sprintf('%s[%s]', newModel.metNames{i}, newModel.comps{newModel.metComps(i)});
end
[~, order] = sort(labels);
metabolites = cell(1, numel(order));
for k = 1:numel(order)
    i = order(k);
    formula = '';
    if isfield(newModel, 'metFormulas'), formula = newModel.metFormulas{i}; end
    charge = NaN;
    if isfield(newModel, 'metCharges'), charge = double(newModel.metCharges(i)); end
    metabolites{k} = struct('metabolite', labels{i}, 'formula', formula, 'charge', charge);
end

nR = numel(newModel.rxns);
reactions = cell(1, nR);
names = cell(1, nR);
for j = 1:nR
    rows = find(newModel.S(:, j));
    [~, o] = sort(labels(rows));
    rows = rows(o);
    entries = cell(1, numel(rows));
    for q = 1:numel(rows)
        entries{q} = struct('metabolite', labels{rows(q)}, 'coefficient', full(double(newModel.S(rows(q), j))));
    end
    clauses = {};
    if ~isempty(newModel.grRules{j})
        dnf = grRuleToDNF(newModel.grRules{j});
        clauses = cell(1, numel(dnf));
        for q = 1:numel(dnf)
            clauses{q} = sort(dnf{q}(:))';
        end
    end
    reactions{j} = struct('name', newModel.rxnNames{j}, 'lower_bound', double(newModel.lb(j)), ...
        'upper_bound', double(newModel.ub(j)), 'clauses', {clauses}, 'stoichiometry', {entries});
    names{j} = newModel.rxnNames{j};
end
[~, o] = sort(names);
reactions = reactions(o);

shortNames = repmat({''}, numel(newModel.genes), 1);
if isfield(newModel, 'geneShortNames'), shortNames = newModel.geneShortNames; end
[~, o] = sort(newModel.genes);
genes = cell(1, numel(o));
for k = 1:numel(o)
    genes{k} = struct('gene', newModel.genes{o(k)}, 'short_name', shortNames{o(k)});
end

results = struct('error', false, 'metabolites', {metabolites}, 'reactions', {reactions}, 'genes', {genes}, ...
    'n_reactions', nR, 'n_metabolites', numel(newModel.mets));
end


function out = as_list(value)
% jsondecode gives a struct array when every element has the same fields and a cell otherwise.
if iscell(value)
    out = value(:)';
else
    out = num2cell(value(:)');
end
end
