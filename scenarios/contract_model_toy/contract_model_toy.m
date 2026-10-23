function results = contract_model_toy(ctx)
% MATLAB side of the contractModel scenario on hand-written models.

inputs = ctx.inputs;
mets = as_list(inputs.mets);
cases = as_list(inputs.cases);
results = struct();
for k = 1:numel(cases)
    c = cases{k};
    model = toy_model(mets, as_list(c.reactions));
    args = {};
    if ~isempty(c.mets)
        args = {'mets', reshape(cellfun(@char, as_list(c.mets), 'UniformOutput', false), [], 1)};
    end
    state = warning('off', 'all');
    [out, removed] = contractModel(model, args{:});
    warning(state);
    [sortedRxns, order] = sort(out.rxns);
    records = cell(1, numel(sortedRxns));
    for r = 1:numel(sortedRxns)
        i = order(r);
        records{r} = struct('id', sortedRxns{r}, 'lb', out.lb(i), 'ub', out.ub(i), 'c', out.c(i), ...
            'rule', {clauses(out.grRules{i})});
    end
    results.(char(c.name)) = struct('reactions', {records}, 'removed', {sort(removed(:))'});
end
end


function out = clauses(rule)
% Top-level OR clauses of a gene rule, without redundant outer parentheses, sorted.
tokens = strsplit(strtrim(strrep(strrep(rule, '(', ' ( '), ')', ' ) ')));
tokens = tokens(~cellfun(@isempty, tokens));
parts = {};
depth = 0;
current = '';
for t = 1:numel(tokens)
    token = tokens{t};
    if strcmp(token, '(')
        depth = depth + 1;
    elseif strcmp(token, ')')
        depth = depth - 1;
    end
    if strcmp(token, 'or') && depth == 0
        parts{end+1} = current; %#ok<AGROW>
        current = '';
    else
        current = [current token]; %#ok<AGROW>
    end
end
if ~isempty(current)
    parts{end+1} = current;
end
for p = 1:numel(parts)
    cl = parts{p};
    while numel(cl) > 1 && cl(1) == '(' && cl(end) == ')' && sum(cl == '(') == 1
        cl = cl(2:end-1);
    end
    parts{p} = cl;
end
out = unique(parts);
if isempty(out)
    out = {};
end
out = reshape(out, 1, []);
end


function model = toy_model(mets, reactions)
n = numel(reactions);
m = numel(mets);
model = struct();
model.id = 'contract_toy';
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
model.lb = zeros(n, 1);
model.ub = zeros(n, 1);
model.c = zeros(n, 1);
model.grRules = cell(n, 1);
S = zeros(m, n);
genes = {};
for j = 1:n
    r = reactions{j};
    model.rxns{j} = char(r.id);
    model.rxnNames{j} = char(r.id);
    model.lb(j) = double(r.lb);
    model.ub(j) = double(r.ub);
    model.c(j) = double(r.c);
    model.grRules{j} = char(r.rule);
    tokens = regexp(model.grRules{j}, '[A-Za-z0-9_]+', 'match');
    genes = [genes, tokens(~ismember(tokens, {'and', 'or'}))]; %#ok<AGROW>
    for e = as_list(r.stoichiometry)
        S(strcmp(model.mets, char(e{1}.met)), j) = double(e{1}.coefficient);
    end
end
model.S = sparse(S);
model.rev = double(model.lb < 0);
model.b = zeros(m, 1);
model.genes = unique(genes(:));
model.rxnGeneMat = sparse(n, numel(model.genes));
for j = 1:n
    tokens = regexp(model.grRules{j}, '[A-Za-z0-9_]+', 'match');
    tokens = tokens(~ismember(tokens, {'and', 'or'}));
    model.rxnGeneMat(j, ismember(model.genes, tokens)) = 1;
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
