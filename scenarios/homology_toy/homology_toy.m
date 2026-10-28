function results = homology_toy(ctx)
% MATLAB side of the getModelFromHomology scenario on one-reaction templates.
%
% Templates convert a_c to b_c; only the gene rules differ. The orthologs are fixed, so the
% model-building logic is the only variable.

cases = as_list(ctx.inputs.cases);
results = struct();
for k = 1:numel(cases)
    c = cases{k};
    specs = as_list(c.templates);
    templates = cell(1, numel(specs));
    blast = [];
    pairs = as_list(c.pairs);
    for t = 1:numel(specs)
        templates{t} = template_model(specs{t});
        list = cell(0, 2);
        for p = 1:numel(pairs)
            if strcmp(char(pairs{p}.template), templates{t}.id)
                list(end+1, :) = {char(pairs{p}.source), char(pairs{p}.target)}; %#ok<AGROW>
            end
        end
        % makeFakeBlastStructure wants at least 10 pairs; the padding pairs name genes that no template has.
        for pad = size(list, 1)+1:10
            list(pad, :) = {sprintf('pad_%s_%d', templates{t}.id, pad), sprintf('npad_%s_%d', templates{t}.id, pad)};
        end
        fake = makeFakeBlastStructure(list, templates{t}.id, char(c.target));
        if isempty(blast)
            blast = fake;
        else
            blast = [blast, fake]; %#ok<AGROW>
        end
    end
    args = {'complexPolicy', char(c.complex_policy)};
    order = as_cellstr(c.preferred_order);
    if ~isempty(order)
        args = [args, {'preferredOrder', order}]; %#ok<AGROW>
    end
    try
        evalc('draft = getModelFromHomology(templates, blast, char(c.target), args{:});');
    catch err
        results.(char(c.name)) = struct('error', err.message);
        continue
    end
    [sortedRxns, idx] = sort(draft.rxns);
    records = cell(1, numel(sortedRxns));
    for r = 1:numel(sortedRxns)
        rule = draft.grRules{idx(r)};
        clauses = {};
        if ~isempty(rule)
            dnf = grRuleToDNF(rule);
            clauses = cell(1, numel(dnf));
            for q = 1:numel(dnf)
                genes = dnf{q};
                clauses{q} = sort(genes(:))';
            end
            keys = cellfun(@(x) strjoin(x, char(0)), clauses, 'UniformOutput', false);
            [~, o] = sort(keys);
            clauses = clauses(o);
        end
        records{r} = struct('reaction', sortedRxns{r}, 'clauses', {clauses});
    end
    results.(char(c.name)) = struct('n_reactions', numel(sortedRxns), 'reactions', {records});
end
end


function model = template_model(spec)
reactions = as_list(spec.reactions);
n = numel(reactions);
model = struct();
model.id = char(spec.id);
model.rxns = cell(n, 1);
model.grRules = cell(n, 1);
model.mets = {'a_c'; 'b_c'};
model.metNames = {'a'; 'b'};
model.S = sparse(repmat([-1; 1], 1, n));
model.lb = zeros(n, 1);
model.ub = 1000 * ones(n, 1);
model.rev = zeros(n, 1);
model.c = zeros(n, 1);
model.b = zeros(2, 1);
model.comps = {'c'};
model.compNames = {'c'};
model.metComps = [1; 1];
genes = {};
for i = 1:n
    model.rxns{i} = char(reactions{i}.id);
    model.grRules{i} = char(reactions{i}.rule);
    tokens = regexp(model.grRules{i}, '[A-Za-z0-9_.]+', 'match');
    tokens = tokens(~ismember(lower(tokens), {'and', 'or'}));
    genes = [genes, tokens]; %#ok<AGROW>
end
model.genes = unique(genes, 'stable')';
model.rxnGeneMat = sparse(n, numel(model.genes));
for i = 1:n
    tokens = regexp(model.grRules{i}, '[A-Za-z0-9_.]+', 'match');
    tokens = tokens(~ismember(lower(tokens), {'and', 'or'}));
    [~, col] = ismember(tokens, model.genes);
    model.rxnGeneMat(i, col) = 1;
end
model.rxnNames = model.rxns;
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
