function results = edge_model_toy(ctx)
% MATLAB side of the edge-case model scenario.
%
% The model is built from scenario.yml, written to YAML and read back, split into irreversible
% reactions, and one gene is removed with the default arguments.

inputs = ctx.inputs;

workdir = tempname;
mkdir(workdir);
cleanup = onCleanup(@() rmdir(workdir, 's'));
path = fullfile(workdir, 'edge.yml');
try
    model = build(inputs);
    evalc('writeYAMLmodel(model, path);');
    lines = splitlines(fileread(path));
    if ~isempty(lines) && isempty(lines{end})
        lines(end) = [];
    end
    keep = ~startsWith(strtrim(lines), {'date:', '- date:'});
    lines = lines(keep);
    evalc('reread = readYAMLmodel(path);');
    rereadSummary = summary(reread);
catch err
    lines = {};
    rereadSummary = struct('error', err.message);
end

try
    evalc('writeYAMLmodel(build(inputs, true), fullfile(workdir, ''dup.yml''));');
    duplicateWritten = true;
catch
    duplicateWritten = false;
end

irrev = convertToIrrev(build(inputs));
evalc('removed = removeGenes(build(inputs), {char(inputs.removed_gene)});');

results.duplicate_names_written = duplicateWritten;
results.yaml_lines = reshape(lines, 1, []);
results.reread = rereadSummary;
results.irreversible = summary(irrev);
results.removed_gene = summary(removed);
end


function model = build(inputs, duplicate)
if nargin < 2
    duplicate = false;
end
reactions = as_list(inputs.reactions);
n = numel(reactions);
names = {};
for i = 1:n
    names = [names, as_cellstr(reactions{i}.substrates), as_cellstr(reactions{i}.products)]; %#ok<AGROW>
end
names = [unique(names, 'stable'), {char(inputs.orphan_metabolite)}];
S = zeros(numel(names), n);
model = struct();
model.id = 'edge';
model.name = 'edge';
model.rxns = cell(n, 1);
model.rxnNames = cell(n, 1);
model.lb = zeros(n, 1);
model.ub = zeros(n, 1);
model.c = zeros(n, 1);
model.grRules = cell(n, 1);
for i = 1:n
    r = reactions{i};
    model.rxns{i} = char(r.id);
    model.rxnNames{i} = char(r.name);
    model.lb(i) = double(r.lb);
    model.ub(i) = double(r.ub);
    model.c(i) = double(r.objective);
    model.grRules{i} = char(r.rule);
    for m = as_cellstr(r.substrates)
        S(strcmp(names, m{1}), i) = -1;
    end
    for m = as_cellstr(r.products)
        S(strcmp(names, m{1}), i) = 1;
    end
end
model.mets = strcat(names(:), '_c');
display = names(:);
if duplicate
    display(ismember(display, {'X1', 'X2'})) = {char(inputs.duplicate_name)};
end
model.metNames = display;
model.metFormulas = repmat({'C'}, numel(names), 1);
model.S = sparse(S);
model.rev = double(model.lb < 0);
model.b = zeros(numel(names), 1);
model.comps = {'c'};
model.compNames = {'cytosol'};
model.metComps = ones(numel(names), 1);
genes = {};
for i = 1:n
    tokens = regexp(model.grRules{i}, '[A-Za-z0-9_.]+', 'match');
    genes = [genes, tokens(~ismember(lower(tokens), {'and', 'or'}))]; %#ok<AGROW>
end
genes = [unique(genes, 'stable'), {char(inputs.orphan_gene)}];
model.genes = genes(:);
model.rxnGeneMat = sparse(n, numel(model.genes));
for i = 1:n
    tokens = regexp(model.grRules{i}, '[A-Za-z0-9_.]+', 'match');
    tokens = tokens(~ismember(lower(tokens), {'and', 'or'}));
    [~, col] = ismember(tokens, model.genes);
    model.rxnGeneMat(i, col) = 1;
end
end


function out = summary(model)
[sortedRxns, order] = sort(model.rxns);
records = cell(1, numel(sortedRxns));
for k = 1:numel(sortedRxns)
    i = order(k);
    clauses = {};
    if ~isempty(model.grRules{i})
        dnf = grRuleToDNF(model.grRules{i});
        clauses = cell(1, numel(dnf));
        for q = 1:numel(dnf)
            clauses{q} = sort(dnf{q}(:))';
        end
        keys = cellfun(@(x) strjoin(x, char(0)), clauses, 'UniformOutput', false);
        [~, o] = sort(keys);
        clauses = clauses(o);
    end
    records{k} = struct('reaction', sortedRxns{k}, 'lower_bound', double(model.lb(i)), ...
        'upper_bound', double(model.ub(i)), 'clauses', {clauses}, 'objective', double(model.c(i)));
end
out.reactions = records;
out.n_metabolites = numel(model.mets);
out.genes = reshape(sort(model.genes), 1, []);
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


function out = as_list(value)
% jsondecode gives a struct array when every element has the same fields and a cell otherwise.
if iscell(value)
    out = value(:)';
else
    out = num2cell(value(:)');
end
end
