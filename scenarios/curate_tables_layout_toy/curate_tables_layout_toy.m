function results = curate_tables_layout_toy(ctx)
% MATLAB side of the curateModelFromTables table layout scenario.
%
% The model has two metabolites, x and y in compartment c, and the reactions R0 (x to y) and R0b (y to x). Each case
% curates it from a set of table files and records the metabolites (name, compartment, formula, charge) and the
% reactions (name, bounds and stoichiometry by metabolite name) so that identifiers, which each side numbers its own
% way, are not compared. The yaml_round_trip cases write the curated model with writeYAMLmodel and read it
% back with readYAMLmodel before recording it.

inputs = ctx.inputs;
cases = as_list(inputs.cases);
results = struct();
for k = 1:numel(cases)
    results.(char(cases{k}.name)) = curate(inputs, cases{k});
end
% One new reaction per gene rule, written to a table of its own; the model holds every gene the rules name.
rules = as_list(inputs.gene_rules);
scratch = tempname;
mkdir(scratch);
cleanup = onCleanup(@() rmdir(scratch, 's'));
for n = 1:numel(rules)
    rxnsFile = fullfile(scratch, sprintf('rxns_%d.tsv', n));
    fid = fopen(rxnsFile, 'w');
    fprintf(fid, 'rxnIdx\trxnNames\tgrRules\n1\tfresh\t%s\n', char(rules{n}));
    fclose(fid);
    c = struct('name', sprintf('gene_rule_%02d', n), 'mets', 'none', 'rxns', rxnsFile, 'coeffs', 'coeffs_fresh.tsv');
    results.(c.name) = curate(inputs, c);
end
rules = as_list(inputs.yaml_round_trip_rules);
for n = 1:numel(rules)
    rxnsFile = fullfile(scratch, sprintf('rxns_yaml_%d.tsv', n));
    fid = fopen(rxnsFile, 'w');
    fprintf(fid, 'rxnIdx\trxnNames\tgrRules\n1\tfresh\t%s\n', char(rules{n}));
    fclose(fid);
    c = struct('name', sprintf('yaml_round_trip_%02d', n), 'mets', 'none', 'rxns', rxnsFile, ...
        'coeffs', 'coeffs_fresh.tsv', 'yaml_file', fullfile(scratch, sprintf('model_%d.yml', n)));
    results.(c.name) = curate(inputs, c);
end
end


function results = curate(inputs, c)
model = struct();
model.id = 'toy';
model.name = 'toy';
model.rxns = {'R0'; 'R0b'};
model.rxnNames = {'existing'; 'existing2'};
model.mets = {'x_c'; 'y_c'};
model.metNames = {'x'; 'y'};
model.metFormulas = {'CH4'; 'C2H6'};
model.metCharges = [1; 2];
model.S = sparse([-1 1; 1 -1]);
model.lb = [0; 0];
model.ub = [1000; 1000];
model.rev = [0; 0];
model.c = [0; 0];
model.b = zeros(2, 1);
model.comps = {'c'};
model.compNames = {'c'};
model.metComps = [1; 1];
modelGenes = as_list(inputs.model_genes);
model.genes = cellfun(@char, modelGenes(:), 'UniformOutput', false);
model.grRules = {''; ''};
model.rxnGeneMat = sparse(2, numel(modelGenes));

folder = char(inputs.tables_dir);
args = {};
names = {'rxns', 'rxnsInfo'; 'coeffs', 'rxnsCoeffs'; 'genes', 'genesInfo'};
for q = 1:size(names, 1)
    if isfield(c, names{q, 1}) && ~strcmp(char(c.(names{q, 1})), 'none')
        file = char(c.(names{q, 1}));
        if ~(numel(file) > 0 && (file(1) == '/' || (numel(file) > 1 && file(2) == ':')))
            file = fullfile(folder, file);
        end
        args = [args, {names{q, 2}, file}]; %#ok<AGROW>
    end
end
metsInfo = 'none';
if isfield(c, 'mets') && ~strcmp(char(c.mets), 'none')
    metsInfo = fullfile(folder, char(c.mets));
end
args = [args, {'metPrefix', char(inputs.met_prefix), 'rxnPrefix', char(inputs.rxn_prefix)}];

try
    evalc('newModel = curateModelFromTables(model, metsInfo, args{:});');
catch err
    kind = error_class(err.message);
    if isfield(c, 'ignore_error_class') && c.ignore_error_class
        kind = '';
    end
    results = struct('error', true, 'error_class', kind);
    return
end
if isfield(c, 'yaml_file')
    evalc('writeYAMLmodel(newModel, c.yaml_file);');
    evalc('newModel = readYAMLmodel(c.yaml_file);');
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
rxnNames = cell(1, nR);
for j = 1:nR
    rows = find(newModel.S(:, j));
    [~, o] = sort(labels(rows));
    rows = rows(o);
    stoich = cell(1, numel(rows));
    for q = 1:numel(rows)
        stoich{q} = struct('metabolite', labels{rows(q)}, 'coefficient', full(double(newModel.S(rows(q), j))));
    end
    ruleGenes = sort(newModel.genes(newModel.rxnGeneMat(j, :) ~= 0));
    reactions{j} = struct('name', newModel.rxnNames{j}, 'lb', bound(newModel.lb(j)), 'ub', bound(newModel.ub(j)), ...
        'genes', {reshape(ruleGenes, 1, [])}, 'stoichiometry', {stoich});
    rxnNames{j} = newModel.rxnNames{j};
end
[~, o] = sort(rxnNames);
reactions = reactions(o);

results = struct('error', false, 'metabolites', {metabolites}, 'reactions', {reactions}, ...
    'model_genes', {reshape(sort(newModel.genes), 1, [])});
end


function value = bound(x)
% Infinite bounds as text, so that both sides write the same JSON.
if x == Inf
    value = 'Inf';
elseif x == -Inf
    value = '-Inf';
else
    value = double(x);
end
end


function out = as_list(value)
% jsondecode gives a struct array when every element has the same fields and a cell otherwise.
if iscell(value)
    out = value(:)';
else
    out = num2cell(value(:)');
end
end


function kind = error_class(message)
% Kind of table error, from the words both sides use in their messages.
words = {'needs its own index', 'repeated_rxn_index'; 'finite whole number', 'charge'; '0, 1, true or false', 'rev'; 'rxnIdx is empty', 'rxn_index_missing'; ...
    'at least 1', 'rxn_index'; 'no matching row', 'reaction_not_matched'; 'not matched', 'reaction_not_matched'; ...
    'has no compartment', 'compartment'; 'must match model.comps', 'compartment'; ...
    'cells but the header has', 'extra_cell'; 'not a valid gene rule', 'gene_rule'; 'more than once', 'repeated_gene'; ...
    'Not all genes', 'unknown_gene'; 'coefficient', 'coefficient'; ...
    'is not a number', 'not_a_number'};
kind = 'other';
for k = 1:size(words, 1)
    if contains(message, words{k, 1})
        kind = words{k, 2};
        return
    end
end
end
