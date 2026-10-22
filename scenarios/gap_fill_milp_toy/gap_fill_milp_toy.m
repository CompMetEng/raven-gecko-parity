function results = gap_fill_milp_toy(ctx)
% MATLAB side of the growth-floor and LP gap-filling scenario on hand-written models.

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
    ids = as_cellstr(c.objective);
    coeffs = double(c.objective_coeffs(:));
    for i = 1:numel(ids)
        draft.c(strcmp(draft.rxns, ids{i})) = coeffs(i);
    end
    try
        if strcmp(char(c.algorithm), 'fastlp')
            [addedRxns, ~, cannotConnect] = gapFillFastLP(draft, template, 'verbose', false);
            out = struct('added_reactions', {sort_cellstr(addedRxns)}, ...
                'cannot_connect', {sort_cellstr(cannotConnect)});
        else
            [addedRxns, ~, exitFlag] = gapFillMILP(draft, template, 'minGrowth', double(c.min_growth), ...
                'bigM', double(c.big_m), 'verbose', false);
            out = struct('added_reactions', {sort_cellstr(addedRxns)}, 'solved', exitFlag == 1);
        end
    catch err
        out = struct('error', err.message);
    end
    results.(char(c.name)) = out;
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
