function results = generate_new_ids_toy(ctx)
% MATLAB side of the generateNewIds scenario on lists of existing identifiers.

cases = as_list(ctx.inputs.cases);
results = struct();
for k = 1:numel(cases)
    c = cases{k};
    ids = reshape(cellfun(@char, as_list(c.ids), 'UniformOutput', false), [], 1);
    model = struct();
    if strcmp(char(c.type), 'rxns')
        model.rxns = ids;
        model.mets = cell(0, 1);
    else
        model.mets = ids;
        model.rxns = cell(0, 1);
    end
    try
        evalc("newIds = generateNewIds(model, char(c.type), char(c.prefix), 'quantity', double(c.quantity), 'numLength', double(c.num_length));");
    catch
        results.(char(c.name)) = struct('error', true);
        continue
    end
    results.(char(c.name)) = struct('ids', {reshape(newIds, 1, [])});
end
end


function out = as_list(value)
% jsondecode gives a struct array when every element has the same fields and a cell otherwise.
if isempty(value)
    out = {};
elseif iscell(value)
    out = value(:)';
else
    out = num2cell(value(:)');
end
end
