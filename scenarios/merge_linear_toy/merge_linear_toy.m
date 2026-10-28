function results = merge_linear_toy(ctx)
% MATLAB side of the mergeLinear scenario on hand-written models.

inputs = ctx.inputs;
cases = as_list(inputs.cases);
results = struct();
for k = 1:numel(cases)
    c = cases{k};
    model = toy_model('toy', as_list(c.model));
    [reduced, origIds, groupIds, reversedRxns] = mergeLinear(model, {});
    origIds = origIds(:);
    groupIds = double(groupIds(:));
    reversedRxns = logical(reversedRxns(:));
    groups = unique(groupIds(groupIds ~= 0));
    members = cell(numel(groups), 1);
    for g = 1:numel(groups)
        members{g} = sort_cellstr(origIds(groupIds == groups(g)));
    end
    % Groups ordered by their first member, as the Python side does.
    if ~isempty(members)
        first = cellfun(@(x) x{1}, members, 'UniformOutput', false);
        [~, order] = sort(first);
        members = members(order);
    end
    records = cell(1, numel(members));
    for g = 1:numel(members)
        records{g} = struct('members', {members{g}});
    end
    results.(char(c.name)) = struct( ...
        'groups', {records}, ...
        'reactions', {sort_cellstr(reduced.rxns)}, ...
        'reversed_reactions', {sort_cellstr(origIds(reversedRxns))});
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
