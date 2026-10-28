function results = task_list_parsing(ctx)
% MATLAB side of the task-list parsing scenario.
%
% parseTaskList returns a struct array holding each bound-carrying field as a
% cell array of names beside two parallel numeric fields (inputs / LBin / UBin).
% run.py gets triples instead, so both are flattened here into the same list of
% records --- neither side's storage should show through into the comparison.
%
% Nothing is sorted: task order, and the order of entries within a task, are
% what the parser made of a sequence of rows.

results = parse_file(ctx.inputs.task_file);

results.edge_cases = struct();
names = ctx.inputs.edge_cases;
if ischar(names)
    names = {names};
end
for k = 1:numel(names)
    name = char(names{k});
    results.edge_cases.(name) = parse_file(ctx.inputs.(name));
end

end


function out = parse_file(path)
% A file RAVEN refuses is recorded as an error without its message.
try
    taskStruct = parseTaskList(path);
catch
    out = struct('error', true);
    return
end

records = cell(1, numel(taskStruct));
for k = 1:numel(taskStruct)
    records{k} = task_record(taskStruct(k));
end

out.error = false;
out.n_tasks = numel(records);
% Cell array of structs, not a struct array: jsonencode turns a 1x1 struct
% array into a bare object instead of a one-element array.
out.tasks = records;
end


function record = task_record(task)
record = struct( ...
    'id', char_field(task, 'id'), ...
    'description', char_field(task, 'description'), ...
    'should_fail', logical(field_or(task, 'shouldFail', false)), ...
    'print_fluxes', logical(field_or(task, 'printFluxes', false)), ...
    'comments', char_field(task, 'comments'), ...
    'inputs', {bounded(task, 'inputs', 'LBin', 'UBin', 'metabolite')}, ...
    'outputs', {bounded(task, 'outputs', 'LBout', 'UBout', 'metabolite')}, ...
    'equations', {bounded(task, 'equations', 'LBequ', 'UBequ', 'equation')}, ...
    'changed', {bounded(task, 'changed', 'LBrxn', 'UBrxn', 'reaction')});
end


function out = bounded(task, nameField, lbField, ubField, nameKey)
% Every key is always present, so a task with no equations is an empty list
% rather than a missing field.
names = field_or(task, nameField, {});
lb = field_or(task, lbField, []);
ub = field_or(task, ubField, []);

names = names(:)';
out = cell(1, numel(names));
for k = 1:numel(names)
    out{k} = struct( ...
        nameKey, char(names{k}), ...
        'lb', bound_at(lb, k), ...
        'ub', bound_at(ub, k));
end
if isempty(out)
    out = {};
end
end


function value = bound_at(bounds, k)
% A bound that is absent where a name is present would be a real difference,
% so it is reported as a NaN rather than quietly defaulted here. (jsonencode
% writes NaN as null and the Python side canonicalises it to "NaN", so the two
% do not match by accident either.)
if numel(bounds) >= k
    value = double(bounds(k));
else
    value = NaN;
end
end


function value = field_or(s, name, fallback)
if isfield(s, name)
    value = s.(name);
else
    value = fallback;
end
end


function value = char_field(s, name)
if isfield(s, name) && ~isempty(s.(name))
    value = char(s.(name));
else
    value = '';
end
end
