function results = ftinit_task_metabolites_toy(ctx)
% MATLAB side of the task-metabolite scenario for ftINIT.
%
% The metabolite T is used only by the reactions R3 and R5, whose genes have low expression, so the
% extraction removes both. The metabolite U has only the dead-end reaction R4, which the first simplification
% removes. Tasks name T and U as outputs. The result records the metabolites and
% reactions of the extracted model and the tasks the prepared data keeps.

inputs = ctx.inputs;
previousSolver = getpref('RAVEN', 'solver', '');
restoreSolver = onCleanup(@() restore_solver(previousSolver));
setRavenSolver(char(inputs.matlab_solver));

model = struct();
model.id = 'toy';
model.mets = {'A_e'; 'B_c'; 'C_c'; 'T_c'; 'U_c'};
model.metNames = {'A'; 'B'; 'C'; 'T'; 'U'};
model.comps = {'c'; 'e'};
model.compNames = {'c'; 'e'};
model.metComps = [2; 1; 1; 1; 1];
model.rxns = {};
model.S = [];
model.rev = [];
model.genes = {'G2'; 'G3'; 'G5'};
model.grRules = {};
model.rxnGeneMat = [];
toAdd = struct();
toAdd.rxns = {'EXA'; 'R1'; 'R2'; 'R3'; 'R5'; 'EXC'; 'R4'};
toAdd.equations = {'=> A[e]'; 'A[e] => B[c]'; 'B[c] <=> C[c]'; 'B[c] => T[c]'; 'T[c] => C[c]'; 'C[c] =>'; 'B[c] => U[c]'};
toAdd.grRules = {''; ''; 'G2'; 'G3'; 'G5'; ''; ''};
evalc('model = addRxns(model, toAdd, 3);');
model.lb = [0; 0; -1000; 0; 0; 0; 0];
model.ub = 1000 * ones(7, 1);
model.c = zeros(7, 1);
model.rxnNames = model.rxns;
model.b = zeros(5, 1);

tasks = parseTaskList(char(inputs.tasks));
evalc('prepData = prepINITModel(model, tasks, ''extComp'', char(inputs.ext_comp));');

genes = sort(model.genes(:));
levels = zeros(numel(genes), 1);
for i = 1:numel(genes)
    levels(i) = double(inputs.levels.(genes{i}));
end
arrayData = struct('genes', {genes}, 'tissues', {{char(inputs.tissue)}}, 'celltypes', {{''}}, ...
    'levels', levels, 'threshold', double(inputs.threshold));
evalc('res = ftINIT(prepData, char(inputs.tissue), [], [], ''transcrData'', arrayData);');

kept = cell(1, numel(prepData.taskStruct));
for q = 1:numel(prepData.taskStruct)
    kept{q} = char(prepData.taskStruct(q).id);
end
results = struct();
results.metabolites = labels(res);
results.reactions = sort(res.rxns(:))';
results.kept_tasks = kept;
results.reference_metabolites = labels(prepData.refModel);
end


function out = labels(m)
out = cell(numel(m.mets), 1);
for i = 1:numel(out)
    out{i} = sprintf('%s[%s]', m.metNames{i}, m.comps{m.metComps(i)});
end
out = sort(out(:))';
end


function restore_solver(previous)
if ~isempty(previous)
    setRavenSolver(previous);
end
end
