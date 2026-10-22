function results = sort_model_toy(ctx)
% MATLAB side of the sortModel scenario on a hand-written model.
%
% Each case runs sortModel with the given flags on a fresh model and records the stoichiometry of the
% reactions whose result does not depend on the random search, the order of the metabolites, the order of
% the reactions when no search ran, the position of the reaction without a subsystem and the order score of
% every subsystem before and after.

inputs = ctx.inputs;
mets = as_list(inputs.mets);
reactions = as_list(inputs.reactions);
cases = as_list(inputs.cases);
results = struct();
for k = 1:numel(cases)
    results.(char(cases{k}.name)) = run_case(mets, reactions, cases{k});
end
end


function out = run_case(mets, reactions, c)
model = toy_model(mets, reactions);
before = order_scores(model);
args = {'sortReversible', logical(c.sort_reversible), 'sortMetName', logical(c.sort_met_name), ...
    'sortReactionOrder', logical(c.sort_reaction_order)};
if ~isempty(c.seed)
    args = [args, {'seed', double(c.seed)}];
end
state = warning('off', 'all');
restoreWarnings = onCleanup(@() warning(state));
model = sortModel(model, args{:});
searched = logical(c.sort_reaction_order);
coefficients = struct();
for j = 1:numel(model.rxns)
    if searched && model.rev(j) ~= 0
        continue
    end
    rows = find(model.S(:, j));
    entries = cell(1, numel(rows));
    for q = 1:numel(rows)
        entries{q} = struct('met', model.mets{rows(q)}, 'coefficient', full(double(model.S(rows(q), j))));
    end
    [~, order] = sort(cellfun(@(e) e.met, entries, 'UniformOutput', false));
    coefficients.(model.rxns{j}) = entries(order);
end
noSub = struct();
for j = 1:numel(model.rxns)
    if isempty(model.subSystems{j})
        noSub.(model.rxns{j}) = j;
    end
end
out = struct('scores_before', before, 'scores_after', order_scores(model), ...
    'metabolite_order', {model.mets(:)'}, 'reactions_sorted', {sort(model.rxns(:))'}, ...
    'coefficients', coefficients, 'no_subsystem_positions', noSub);
if ~searched
    out.reaction_order = model.rxns(:)';
end
end


function scores = order_scores(model)
% The number of metabolites that are produced before they are consumed, per subsystem with 2 to 250 reactions.
scores = struct();
names = {};
for j = 1:numel(model.rxns)
    names = [names, model.subSystems{j}(:)']; %#ok<AGROW>
end
names = unique(names);
for i = 1:numel(names)
    rxns = find(cellfun(@(s) any(strcmp(s, names{i})), model.subSystems));
    if numel(rxns) < 2 || numel(rxns) > 250
        continue
    end
    subS = model.S(:, rxns);
    subS = subS(any(subS, 2), :);
    s = [subS, ones(size(subS, 1), 1), -ones(size(subS, 1), 1)];
    r1 = arrayfun(@(x) find(s(x, :) > 0, 1, 'first'), 1:size(s, 1));
    r2 = arrayfun(@(x) find(s(x, :) < 0, 1, 'first'), 1:size(s, 1));
    scores.(names{i}) = sum(r1 < r2);
end
end


function model = toy_model(mets, reactions)
n = numel(reactions);
m = numel(mets);
model = struct();
model.id = 'sort_toy';
model.mets = cell(m, 1);
model.metNames = cell(m, 1);
compartments = cell(m, 1);
for i = 1:m
    model.mets{i} = char(mets{i}.id);
    model.metNames{i} = char(mets{i}.name);
    compartments{i} = char(mets{i}.compartment);
end
[model.comps, ~, model.metComps] = unique(compartments, 'stable');
model.comps = model.comps(:)';
model.compNames = model.comps;
model.metComps = model.metComps(:);
model.rxns = cell(n, 1);
model.rxnNames = cell(n, 1);
model.subSystems = cell(n, 1);
model.lb = zeros(n, 1);
model.ub = zeros(n, 1);
S = zeros(m, n);
for j = 1:n
    r = reactions{j};
    model.rxns{j} = char(r.id);
    model.rxnNames{j} = char(r.id);
    if isempty(r.subsystem)
        model.subSystems{j} = {};
    else
        model.subSystems{j} = {char(r.subsystem)};
    end
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
