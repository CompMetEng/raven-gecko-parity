function results = simplify_order_toy(ctx)
% MATLAB side of the simplification-order scenario on hand-written models.
%
% Same shape as run.py: reactions, metabolites and genes in model order.

inputs = ctx.inputs;
previousSolver = getpref('RAVEN', 'solver', '');
restoreSolver = onCleanup(@() restore_solver(previousSolver));
setRavenSolver(char(inputs.matlab_solver));

modeNames = struct('groupLinear', 'groupLinear', 'constrainReversible', 'constrainReversible', ...
    'deleteZeroInterval', 'deleteZeroInterval', 'deleteDuplicates', 'deleteDuplicates', ...
    'deleteInaccessible', 'deleteInaccessible');

cases = as_list(inputs.cases);
results = struct();
for k = 1:numel(cases)
    c = cases{k};
    model = toy_model('toy', as_list(c.reactions), as_cellstr(c.extra_metabolites));
    args = {'deleteUnconstrained', false, 'suppressWarnings', true, 'reservedRxns', as_cellstr(c.reserved)};
    modes = as_cellstr(c.modes);
    for m = 1:numel(modes)
        args = [args, {modeNames.(modes{m}), true}]; %#ok<AGROW>
    end
    reduced = simplifyModel(model, args{:});
    bounds = cell(1, numel(reduced.rxns));
    entries = {};
    for j = 1:numel(reduced.rxns)
        bounds{j} = struct('reaction', reduced.rxns{j}, 'lower_bound', double(reduced.lb(j)), ...
            'upper_bound', double(reduced.ub(j)));
        rows = find(reduced.S(:, j));
        [~, o] = sort(reduced.mets(rows));
        rows = rows(o);
        for q = 1:numel(rows)
            entries{end+1} = struct('reaction', reduced.rxns{j}, 'metabolite', reduced.mets{rows(q)}, ... %#ok<AGROW>
                'coefficient', full(double(reduced.S(rows(q), j))));
        end
    end
    results.(char(c.name)) = struct( ...
        'reactions', {row(reduced.rxns)}, ...
        'metabolites', {row(reduced.mets)}, ...
        'genes', {row(reduced.genes)}, ...
        'bounds', {bounds}, ...
        'stoichiometry', {entries});
end
end


function model = toy_model(id, reactions, extraMetabolites)
n = numel(reactions);
mets = {};
for i = 1:n
    mets = [mets, as_cellstr(reactions{i}.substrates), as_cellstr(reactions{i}.products)]; %#ok<AGROW>
end
mets = [unique(mets, 'stable'), extraMetabolites];
S = zeros(numel(mets), n);
genes = {};
rules = repmat({''}, n, 1);
model = struct();
model.id = id;
model.rxns = cell(n, 1);
model.lb = zeros(n, 1);
model.ub = zeros(n, 1);
for i = 1:n
    r = reactions{i};
    model.rxns{i} = char(r.id);
    model.lb(i) = double(r.lb);
    model.ub(i) = double(r.ub);
    for m = as_cellstr(r.substrates)
        S(strcmp(mets, m{1}), i) = -1;
    end
    for m = as_cellstr(r.products)
        S(strcmp(mets, m{1}), i) = 1;
    end
    if ~isempty(r.rule)
        rules{i} = char(r.rule);
        genes = [genes, {char(r.rule)}]; %#ok<AGROW>
    end
end
genes = unique(genes, 'stable');
rxnGeneMat = zeros(n, numel(genes));
for i = 1:n
    if ~isempty(rules{i})
        rxnGeneMat(i, strcmp(genes, rules{i})) = 1;
    end
end
model.mets = mets(:);
model.metNames = mets(:);
model.S = sparse(S);
model.rev = double(model.lb < 0 & model.ub > 0);
model.c = zeros(n, 1);
model.b = zeros(numel(mets), 1);
model.comps = {'c'};
model.compNames = {'c'};
model.metComps = ones(numel(mets), 1);
model.rxnNames = model.rxns;
model.genes = genes(:);
model.grRules = rules;
model.rxnGeneMat = sparse(rxnGeneMat);
end


function restore_solver(previous)
if ~isempty(previous)
    setRavenSolver(previous);
end
end


function out = row(values)
out = reshape(values, 1, []);
if isempty(out)
    out = {};
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


function out = as_list(value)
% jsondecode gives a struct array when every element has the same fields and a cell otherwise.
if iscell(value)
    out = value(:)';
else
    out = num2cell(value(:)');
end
end
