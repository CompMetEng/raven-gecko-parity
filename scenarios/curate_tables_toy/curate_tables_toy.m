function results = curate_tables_toy(ctx)
% MATLAB side of the curateModelFromTables scenario.
%
% The model has two empty reactions; the tables give them metabolites, rules and annotations.
% Content is recorded by metabolite name and compartment, identifiers separately.

inputs = ctx.inputs;
cases = as_list(inputs.cases);
results = struct();
for k = 1:numel(cases)
    results.(char(cases{k}.name)) = curate(inputs, as_cellstr(cases{k}.reactions));
end
end


function results = curate(inputs, empty)
n = numel(empty);
model = struct();
model.id = 'toy';
model.rxns = empty(:);
model.rxnNames = repmat({''}, n, 1);
model.mets = cell(0, 1);
model.metNames = cell(0, 1);
model.S = sparse(0, n);
model.lb = zeros(n, 1);
model.ub = 1000 * ones(n, 1);
model.rev = zeros(n, 1);
model.c = zeros(n, 1);
model.b = zeros(0, 1);
model.comps = {'c'};
model.compNames = {'c'};
model.metComps = zeros(0, 1);
model.genes = cell(0, 1);
model.grRules = repmat({''}, n, 1);
model.rxnGeneMat = sparse(n, 0);

try
    evalc(['newModel = curateModelFromTables(model, char(inputs.mets_table), ' ...
        '''genesInfo'', char(inputs.genes_table), ''rxnsInfo'', char(inputs.rxns_table), ' ...
        '''rxnsCoeffs'', char(inputs.coeffs_table), ''metPrefix'', char(inputs.met_prefix), ' ...
        '''rxnPrefix'', char(inputs.rxn_prefix));']);
catch err
    results = struct('error', err.message);
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
keys = cell(1, nR);
for j = 1:nR
    rows = find(newModel.S(:, j));
    entries = cell(1, numel(rows));
    [~, o] = sort(labels(rows));
    rows = rows(o);
    for q = 1:numel(rows)
        entries{q} = struct('metabolite', labels{rows(q)}, 'coefficient', full(double(newModel.S(rows(q), j))));
    end
    clauses = {};
    if ~isempty(newModel.grRules{j})
        dnf = grRuleToDNF(newModel.grRules{j});
        clauses = cell(1, numel(dnf));
        for c = 1:numel(dnf)
            clauses{c} = sort(dnf{c}(:))';
        end
    end
    reactions{j} = struct('name', newModel.rxnNames{j}, 'lower_bound', double(newModel.lb(j)), ...
        'upper_bound', double(newModel.ub(j)), 'clauses', {clauses}, 'stoichiometry', {entries});
    clauseText = cellfun(@(c) strjoin(c, ','), clauses, 'UniformOutput', false);
    keys{j} = sprintf('%s|%s|%s', newModel.rxnNames{j}, strjoin(labels(rows)', ','), strjoin(clauseText, ';'));
end
[~, o] = sort(keys);
reactions = reactions(o);

shortNames = repmat({''}, numel(newModel.genes), 1);
if isfield(newModel, 'geneShortNames'), shortNames = newModel.geneShortNames; end
[~, o] = sort(newModel.genes);
genes = cell(1, numel(o));
for k = 1:numel(o)
    genes{k} = struct('gene', newModel.genes{o(k)}, 'short_name', shortNames{o(k)});
end

results = struct();
results.metabolites = metabolites;
results.reactions = reactions;
results.genes = genes;
results.n_reactions = nR;
results.n_metabolites = numel(newModel.mets);
results.reaction_ids = reshape(sort(newModel.rxns), 1, []);
results.metabolite_ids = reshape(sort(newModel.mets), 1, []);
end


function out = as_list(value)
% jsondecode gives a struct array when every element has the same fields and a cell otherwise.
if iscell(value)
    out = value(:)';
else
    out = num2cell(value(:)');
end
end


function out = as_cellstr(value)
if isempty(value)
    out = {};
elseif ischar(value)
    out = {value};
else
    out = reshape(cellfun(@char, value, 'UniformOutput', false), 1, []);
end
end
