function results = curate_tables_annotation_toy(ctx)
% MATLAB side of the curateModelFromTables annotation scenario.
%
% The model has two metabolites with a KEGG compound annotation, two genes with a short name and an NCBI gene
% annotation, and two reactions with a KEGG reaction annotation and a confidence score. Each case curates it from
% a set of tables and records annotations, charges and scores by name, so that identifiers, which each side
% numbers its own way, are not compared.

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
model.rxns = {'R0'; 'R0b'};
model.rxnNames = {'existing'; 'existing2'};
model.mets = {'x_c'; 'y_c'};
model.metNames = {'x'; 'y'};
model.metFormulas = {'CH4'; 'C2H6'};
model.metCharges = [1; 2];
model.metMiriams = {struct('name', {{'kegg.compound'}}, 'value', {{'C1'}}); struct('name', {{'kegg.compound'}}, 'value', {{'C2'}})};
model.S = sparse([-1 1; 1 -1]);
model.lb = [0; 0];
model.ub = [1000; 1000];
model.rev = [0; 0];
model.c = [0; 0];
model.b = zeros(2, 1);
model.comps = {'c'};
model.compNames = {'c'};
model.metComps = [1; 1];
model.genes = {'G1'; 'G2'};
model.geneShortNames = {'one'; 'two'};
model.geneMiriams = {struct('name', {{'ncbigene'}}, 'value', {{'11'}}); struct('name', {{'ncbigene'}}, 'value', {{'22'}})};
model.grRules = {'G1 or G2'; ''};
model.rxnGeneMat = sparse([1 1; 0 0]);
model.rxnMiriams = {struct('name', {{'kegg.reaction'}}, 'value', {{'R1'}}); struct('name', {{'kegg.reaction'}}, 'value', {{'R2'}})};
model.rxnConfidenceScores = [2; 4];

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
catch err
    results = struct('error', true, 'error_class', error_class(err.message));
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
    metabolites{k} = struct('metabolite', labels{i}, 'formula', formula, 'charge', charge, ...
        'kegg_compound', annotation(newModel, 'metMiriams', i, 'kegg.compound'));
end

nR = numel(newModel.rxns);
reactions = cell(1, nR);
names = cell(1, nR);
for j = 1:nR
    score = NaN;
    if isfield(newModel, 'rxnConfidenceScores'), score = double(newModel.rxnConfidenceScores(j)); end
    reactions{j} = struct('name', newModel.rxnNames{j}, ...
        'kegg_reaction', annotation(newModel, 'rxnMiriams', j, 'kegg.reaction'), ...
        'metacyc', annotation(newModel, 'rxnMiriams', j, 'metacyc'), 'confidence_score', score, ...
        'lb', bound(newModel.lb(j)), 'ub', bound(newModel.ub(j)));
    names{j} = newModel.rxnNames{j};
end
[~, o] = sort(names);
reactions = reactions(o);

shortNames = repmat({''}, numel(newModel.genes), 1);
if isfield(newModel, 'geneShortNames'), shortNames = newModel.geneShortNames; end
[~, o] = sort(newModel.genes);
genes = cell(1, numel(o));
for k = 1:numel(o)
    genes{k} = struct('gene', newModel.genes{o(k)}, 'short_name', shortNames{o(k)}, ...
        'ncbigene', annotation(newModel, 'geneMiriams', o(k), 'ncbigene'));
end

results = struct('error', false, 'metabolites', {metabolites}, 'reactions', {reactions}, 'genes', {genes});
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


function value = annotation(model, field, index, namespace)
% The value of one MIRIAM namespace for one entity, empty when the entity has none.
value = '';
if ~isfield(model, field) || index > numel(model.(field)) || isempty(model.(field){index})
    return
end
miriam = model.(field){index};
hit = strcmp(miriam.name, namespace);
if any(hit)
    value = strjoin(miriam.value(hit), '; ');
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
words = {'finite whole number', 'charge'; '0, 1, true or false', 'rev'; 'at least 1', 'rxn_index'; ...
    'coefficient', 'coefficient'; 'is not a number', 'not_a_number'};
kind = 'other';
for k = 1:size(words, 1)
    if contains(message, words{k, 1})
        kind = words{k, 2};
        return
    end
end
end
