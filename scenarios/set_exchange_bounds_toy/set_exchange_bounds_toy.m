function results = set_exchange_bounds_toy(ctx)
% MATLAB side of the setExchangeBounds scenario on a hand-written model.
%
% Each case applies setExchangeBounds to a fresh copy of the model and records the bounds of every reaction and
% the metabolites that were not used.

inputs = ctx.inputs;
cases = as_list(inputs.cases);
results = struct();
for k = 1:numel(cases)
    c = cases{k};
    model = toy_model(inputs, as_list(c.reactions));
    args = {'closeOthers', logical(c.close_others), 'mediaOnly', logical(c.media_only)};
    if ~isempty(c.mets)
        args = [args, {'mets', reshape(cellfun(@char, as_list(c.mets), 'UniformOutput', false), [], 1)}]; %#ok<AGROW>
    end
    if ~isempty(c.lb)
        args = [args, {'lb', double(c.lb)}]; %#ok<AGROW>
    end
    if ~isempty(c.ub)
        args = [args, {'ub', double(c.ub)}]; %#ok<AGROW>
    end
    state = warning('off', 'all');
    try
        [out, unused] = setExchangeBounds(model, args{:});
    catch
        warning(state);
        results.(char(c.name)) = struct('error', true);
        continue
    end
    warning(state);
    [sortedRxns, order] = sort(out.rxns);
    bounds = cell(1, numel(sortedRxns));
    for r = 1:numel(sortedRxns)
        bounds{r} = struct('reaction', sortedRxns{r}, 'lb', out.lb(order(r)), 'ub', out.ub(order(r)));
    end
    unused = sort(unused(:))';
    results.(char(c.name)) = struct('bounds', {bounds}, 'unused', {unused});
end
end


function model = toy_model(inputs, reactions)
mets = as_list(inputs.mets);
n = numel(reactions);
m = numel(mets);
model = struct();
model.id = 'exchange_toy';
model.mets = cell(m, 1);
model.metNames = cell(m, 1);
compartments = cell(m, 1);
for i = 1:m
    model.mets{i} = char(mets{i}.id);
    model.metNames{i} = char(mets{i}.name);
    compartments{i} = char(mets{i}.compartment);
end
compIds = fieldnames(inputs.compartments)';
model.comps = compIds;
model.compNames = cellfun(@(id) char(inputs.compartments.(id)), compIds, 'UniformOutput', false);
[~, model.metComps] = ismember(compartments, compIds);
model.rxns = cell(n, 1);
model.rxnNames = cell(n, 1);
model.lb = zeros(n, 1);
model.ub = zeros(n, 1);
S = zeros(m, n);
for j = 1:n
    r = reactions{j};
    model.rxns{j} = char(r.id);
    model.rxnNames{j} = char(r.id);
    model.lb(j) = double(r.lb);
    model.ub(j) = double(r.ub);
    for e = as_list(r.stoichiometry)
        S(strcmp(model.mets, char(e{1}.met)), j) = double(e{1}.coefficient);
    end
end
model.S = sparse(S);
model.rev = double(model.lb < 0);
model.c = zeros(n, 1);
model.b = zeros(m, 1);
model.genes = cell(0, 1);
model.grRules = repmat({''}, n, 1);
model.rxnGeneMat = sparse(n, 0);
end


function out = as_list(value)
% jsondecode gives a struct array when every element has the same fields and a cell otherwise.
if iscell(value)
    out = value(:)';
else
    out = num2cell(value(:)');
end
end
