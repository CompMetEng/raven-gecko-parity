function results = init_task_metabolites_toy(ctx)
% MATLAB side of the task-metabolite scenario for prepINITModel.
%
% A model has the dead-end reaction R3 (B to T) as the only reaction of T, a second metabolite T2 with the same
% name and compartment as T whose only reaction is the dead end R4, and tasks that name T as an output (MAKE_T)
% and in another letter case (MAKE_T_LOWER). The result records the metabolites and reactions of the prepared reference model, the tasks the
% prepared data keeps and the essential reactions.

inputs = ctx.inputs;
model = struct();
model.id = 'toy';
model.rxns = {'EXA'; 'R1'; 'R2'; 'R3'; 'EXC'; 'R4'};
model.rxnNames = model.rxns;
model.mets = {'A_e'; 'B_c'; 'C_c'; 'T_c'; 'T2_c'};
model.metNames = {'A'; 'B'; 'C'; 'T'; 'T'};
model.comps = {'c'; 'e'};
model.compNames = {'c'; 'e'};
model.metComps = [2; 1; 1; 1; 1];
model.S = sparse([ 1 -1  0  0  0  0; ...
                   0  1 -1 -1  0 -1; ...
                   0  0  1  0 -1  0; ...
                   0  0  0  1  0  0; ...
                   0  0  0  0  0  1]);
model.lb = [0; 0; -1000; 0; 0; 0];
model.ub = 1000 * ones(6, 1);
model.rev = [0; 0; 1; 0; 0; 0];
model.c = zeros(6, 1);
model.b = zeros(5, 1);
model.genes = cell(0, 1);
model.grRules = repmat({''}, 6, 1);
model.rxnGeneMat = sparse(6, 0);

tasks = parseTaskList(char(inputs.tasks));
evalc('prepData = prepINITModel(model, tasks, ''extComp'', char(inputs.ext_comp));');
ref = prepData.refModel;

labels = cell(numel(ref.mets), 1);
for i = 1:numel(labels)
    labels{i} = sprintf('%s[%s]', ref.metNames{i}, ref.comps{ref.metComps(i)});
end
kept = cell(1, numel(prepData.taskStruct));
for q = 1:numel(prepData.taskStruct)
    kept{q} = char(prepData.taskStruct(q).id);
end

results = struct();
results.metabolites = sort(labels(:))';
results.reactions = sort(ref.rxns(:))';
results.kept_tasks = kept;
results.essential_reactions = sort(prepData.essentialRxns(:))';
end
