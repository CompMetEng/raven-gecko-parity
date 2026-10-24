function results = init_output_model_smallyeast(ctx)
% MATLAB side of the prepINITModel and ftINIT scenario on smallYeast with a task list.
%
% The template is read, its boundary reactions are opened and a task list is
% prepared. ftINIT then extracts a model from gene levels that alternate between a score of +1 and -1 in
% sorted gene order. The result records the tasks the prepared data keeps, the merged essential reactions,
% the number of reference-model reactions whose bounds or stoichiometry differ from the template, and the
% bounds and equation of three reactions in the extracted model, and the number of genes of the reference
% model and the extracted model. A task list with a repeated task id is refused.

inputs = ctx.inputs;
previousSolver = getpref('RAVEN', 'solver', '');
restoreSolver = onCleanup(@() restore_solver(previousSolver));
setRavenSolver(char(inputs.matlab_solver));

model = readYAMLmodel(char(inputs.model));
model = setParam(model, 'ub', as_cellstr(inputs.open_bounds), double(inputs.open_ub));
templateRxns = model.rxns;
template = cell(numel(templateRxns), 1);
for j = 1:numel(templateRxns)
    template{j} = struct('lb', model.lb(j), 'ub', model.ub(j), 'eq', {equation(model, j)});
end

duplicateError = false;
try
    evalc('prepINITModel(model, parseTaskList(char(inputs.tasks_duplicate)), ''extComp'', char(inputs.ext_comp));');
catch
    duplicateError = true;
end

tasks = parseTaskList(char(inputs.tasks));
evalc('prepData = prepINITModel(model, tasks, ''extComp'', char(inputs.ext_comp));');

genes = sort(model.genes(:));
levels = zeros(numel(genes), 1);
for i = 1:numel(genes)
    if mod(i, 2) == 1
        levels(i) = double(inputs.gene_level_up);
    else
        levels(i) = double(inputs.gene_level_down);
    end
end
arrayData = struct('genes', {genes}, 'tissues', {{char(inputs.tissue)}}, 'celltypes', {{''}}, ...
    'levels', levels, 'threshold', double(inputs.threshold));
evalc('res = ftINIT(prepData, char(inputs.tissue), [], [], ''transcrData'', arrayData);');

ref = prepData.refModel;
nDiff = 0;
for j = 1:numel(ref.rxns)
    k = find(strcmp(templateRxns, ref.rxns{j}));
    t = template{k};
    mine = equation(ref, j);
    sameEq = numel(mine) == numel(t.eq) && isequal(cellfun(@(a) a.met, mine, 'UniformOutput', false), cellfun(@(a) a.met, t.eq, 'UniformOutput', false));
    if sameEq
        sameEq = all(abs(cellfun(@(a) a.coefficient, mine) - cellfun(@(a) a.coefficient, t.eq)) <= 1e-12);
    end
    if ~(sameEq && ref.lb(j) == t.lb && ref.ub(j) == t.ub)
        nDiff = nDiff + 1;
    end
end

wanted = as_cellstr(inputs.reported_reactions);
reported = cell(1, numel(wanted));
for q = 1:numel(wanted)
    j = find(strcmp(res.rxns, wanted{q}));
    if isempty(j)
        reported{q} = struct('reaction', wanted{q}, 'in_model', false, 'lower_bound', 0, 'upper_bound', 0, 'equation', {{}});
    else
        reported{q} = struct('reaction', wanted{q}, 'in_model', true, 'lower_bound', double(res.lb(j)), ...
            'upper_bound', double(res.ub(j)), 'equation', {equation(res, j)});
    end
end

kept = cell(1, numel(prepData.taskStruct));
for q = 1:numel(prepData.taskStruct)
    kept{q} = struct('id', char(prepData.taskStruct(q).id), 'description', char(prepData.taskStruct(q).description));
end

results = struct();
results.n_tasks_parsed = numel(tasks);
results.kept_tasks = kept;
results.essential_reactions = sort_cellstr(prepData.essentialRxns);
results.duplicate_task_ids_refused = duplicateError;
results.n_ref_reactions = numel(ref.rxns);
results.n_ref_genes = numel(ref.genes);
results.n_output_genes = numel(res.genes);
results.n_ref_reactions_differing_from_template = nDiff;
results.n_output_reactions = numel(res.rxns);
results.output_reactions = sort_cellstr(res.rxns);
results.reported = reported;
end


function out = equation(model, j)
rows = find(model.S(:, j));
[~, order] = sort(model.mets(rows));
rows = rows(order);
out = cell(1, numel(rows));
for q = 1:numel(rows)
    out{q} = struct('met', model.mets{rows(q)}, 'coefficient', full(double(model.S(rows(q), j))));
end
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
