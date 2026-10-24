function results = gap_fill_toy(ctx)
% MATLAB side of the connectivity gap-filling scenario on hand-written models.
%
% Defaults are used for everything a case does not name, so a difference in a default shows up here.

inputs = ctx.inputs;
previousSolver = getpref('RAVEN', 'solver', '');
restoreSolver = onCleanup(@() restore_solver(previousSolver));
setRavenSolver(char(inputs.matlab_solver));

cases = as_list(inputs.cases);
results = struct();
for k = 1:numel(cases)
    c = cases{k};
    draft = toy_model('draft', as_list(c.draft));
    template = toy_model('template', as_list(c.template));
    args = {'allowNetProduction', logical(c.net_production), 'supressWarnings', true};
    if ~isempty(c.scores)
        args = [args, {'rxnScores', double(c.scores(:))}]; %#ok<AGROW>
    end
    try
        [newConnected, cannotConnect, addedRxns] = fillGaps(draft, {template}, args{:});
    catch err
        results.(char(c.name)) = struct('error', err.message);
        continue
    end
    results.(char(c.name)) = struct( ...
        'added_reactions', {sort_cellstr(addedRxns)}, ...
        'newly_connected', {sort_cellstr(newConnected)}, ...
        'cannot_connect', {sort_cellstr(cannotConnect)});
end
end


function model = toy_model(id, reactions)
n = numel(reactions);
mets = {};
for i = 1:n
    mets = [mets, as_cellstr(reactions{i}.substrates), as_cellstr(reactions{i}.products)]; %#ok<AGROW>
end
[mets, ~, ~] = unique(mets, 'stable');
S = zeros(numel(mets), n);
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
end
model.mets = mets(:);
model.metNames = mets(:);
model.S = sparse(S);
model.rev = double(model.lb < 0);
model.c = zeros(n, 1);
model.b = zeros(numel(mets), 1);
model.comps = {'c'};
model.compNames = {'c'};
model.metComps = ones(numel(mets), 1);
model.rxnNames = model.rxns;
model.genes = cell(0, 1);
model.grRules = repmat({''}, n, 1);
model.rxnGeneMat = sparse(n, 0);
end


function restore_solver(previous)
if ~isempty(previous)
    setRavenSolver(previous);
end
end


function out = sort_cellstr(values)
if isempty(values)
    out = {};
    return
end
out = reshape(sort(values(:)), 1, []);
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
