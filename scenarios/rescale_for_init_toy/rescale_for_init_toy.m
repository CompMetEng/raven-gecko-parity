function results = rescale_for_init_toy(ctx)
% MATLAB side of the rescaleModelForINIT scenario.

inputs = ctx.inputs;
rxnSpecs = as_list(inputs.reactions);
cases = as_list(inputs.cases);
results = struct();
for k = 1:numel(cases)
    c = cases{k};
    model = struct();
    model.rxns = cell(numel(rxnSpecs), 1);
    model.lb = zeros(numel(rxnSpecs), 1);
    model.ub = zeros(numel(rxnSpecs), 1);
    metIds = {};
    entries = cell(numel(rxnSpecs), 1);
    for r = 1:numel(rxnSpecs)
        s = rxnSpecs{r};
        model.rxns{r} = char(s.id);
        model.lb(r) = double(s.lb);
        model.ub(r) = double(s.ub);
        names = fieldnames(s.coefficients);
        entries{r} = {names, cellfun(@(n) double(s.coefficients.(n)), names)};
        metIds = union(metIds, names);
    end
    metIds = sort(metIds(:));
    model.mets = metIds;
    model.S = zeros(numel(metIds), numel(rxnSpecs));
    for r = 1:numel(rxnSpecs)
        [~, idx] = ismember(entries{r}{1}, metIds);
        model.S(idx, r) = entries{r}{2};
    end
    model.S = sparse(model.S);
    out = rescaleModelForINIT(model, 'maxStoichDiff', double(c.max_stoich_diff));

    records = cell(1, numel(rxnSpecs));
    for r = 1:numel(rxnSpecs)
        idx = find(out.S(:, r));
        coefs = cell(1, numel(idx));
        for j = 1:numel(idx)
            coefs{j} = struct('met', out.mets{idx(j)}, 'coefficient', full(out.S(idx(j), r)));
        end
        records{r} = struct('id', out.rxns{r}, 'coefficients', {coefs}, ...
            'lb', double(out.lb(r)), 'ub', double(out.ub(r)));
    end
    results.(char(c.name)) = records;
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
