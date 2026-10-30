function results = sensitivity_tuning_edge_cases_ectestgem(ctx)
% MATLAB side of the sensitivity-tuning edge cases, mirroring run.py.

previousSolver = getpref('RAVEN', 'solver', '');
restoreSolver = onCleanup(@() restore_solver(previousSolver));
setRavenSolver(char(ctx.inputs.matlab_solver));

adapter = ModelAdapterManager.getAdapter(ctx.inputs.adapter_matlab);
model = loadConventionalGEM('modelAdapter', adapter);

tuneAdapter = adapter;
tuneAdapter.params.bioRxn = char(ctx.inputs.bio_rxn);

inputs = ctx.inputs;
results = struct();

ecModel = build_model(model, adapter, inputs.order_kcats, struct(), false, []);
[results.order, ecModel] = tune(ecModel, tuneAdapter, inputs.order_desired_growth, {});
results.rerun = tune(ecModel, tuneAdapter, inputs.rerun_desired_growth, {});

names = {'whole', 'long'};
for i = 1:numel(names)
    ecModel = build_model(model, adapter, inputs.(['notes_' names{i} '_kcats']), struct(), false, []);
    results.(['notes_' names{i}]) = tune(ecModel, tuneAdapter, inputs.(['notes_' names{i} '_desired_growth']), {});
end

ecModel = build_model(model, adapter, inputs.ceiling_kcats, struct(), false, inputs.ceiling_pool);
results.ceiling = tune(ecModel, tuneAdapter, inputs.ceiling_desired_growth, {});

ecModel = build_model(model, adapter, inputs.ignore_kcats, struct(), false, []);
results.ignore_unknown_id = tune(ecModel, tuneAdapter, inputs.ignore_desired_growth, {char(inputs.ignore_unknown_id)});
ecModel = build_model(model, adapter, inputs.ignore_kcats, struct(), false, []);
results.ignore_lower_case_id = tune(ecModel, tuneAdapter, inputs.ignore_desired_growth, {char(inputs.ignore_lower_case_id)});

ecModel = build_model(model, adapter, inputs.order_kcats, struct('R5', [0 0]), false, []);
results.zero_growth = tune(ecModel, tuneAdapter, inputs.zero_growth_desired_growth, {});

ecModel = build_model(model, adapter, inputs.light_kcats, struct(), true, []);
light = tune(ecModel, tuneAdapter, inputs.light_desired_growth, {});
results.light.raised = light.raised;
if light.raised
    results.light.n_tuned = 0;
else
    results.light.n_tuned = numel(light.tuned);
end

results.truncate_ties = truncate_ties(inputs);
end


function ecModel = build_model(model, adapter, kcats, bounds, light, pool)
rxnIds = fieldnames(bounds);
for i = 1:numel(rxnIds)
    idx = strcmp(model.rxns, rxnIds{i});
    model.lb(idx) = bounds.(rxnIds{i})(1);
    model.ub(idx) = bounds.(rxnIds{i})(2);
end
for id = {'R2', 'R4'}
    idx = strcmp(model.rxns, id{1});
    model.lb(idx) = 0; model.ub(idx) = 0;
end

ecModel = makeEcModel(model, light, adapter);
ecModel = getECfromGEM(ecModel);

kcatRxns = fieldnames(kcats);
for i = 1:numel(ecModel.ec.rxns)
    base = ecModel.ec.rxns{i};
    if light
        base = base(5:end);
    end
    j = find(strcmp(kcatRxns, base));
    if ~isempty(j)
        ecModel.ec.kcat(i) = kcats.(kcatRxns{j});
        ecModel.ec.source{i} = 'manual';
    end
end

ecModel = applyKcatConstraints(ecModel);
if isempty(pool)
    ecModel = setProtPoolSize(ecModel, [], [], [], adapter);
else
    ecModel = setProtPoolSize(ecModel, double(pool.p_tot), double(pool.f), double(pool.sigma), adapter);
end
end


function [out, tunedModel] = tune(ecModel, tuneAdapter, desired, protToIgnore)
tunedModel = ecModel;
try
    [tunedModel, tunedKcats] = sensitivityTuning(ecModel, ...
        'desiredGrowthRate', double(desired), 'modelAdapter', tuneAdapter, ...
        'protToIgnore', protToIgnore, 'verbose', false);
catch
    out = struct('raised', true);
    return
end

n = numel(tunedKcats.rxns);
tuned = cell(1, n);
for i = 1:n
    tuned{i} = struct('reaction', tunedKcats.rxns{i}, 'enzymes', tunedKcats.enzymes{i}, ...
        'old_kcat', num_value(tunedKcats.oldKcat(i)), 'new_kcat', num_value(tunedKcats.newKcat(i)), ...
        'source', tunedKcats.source{i});
end
[~, order] = sort(cellfun(@(t) t.reaction, tuned, 'UniformOutput', false));

sol = solveLP(tunedModel);
[~, ecOrder] = sort(tunedModel.ec.rxns);
ec = cell(1, numel(ecOrder));
for k = 1:numel(ecOrder)
    i = ecOrder(k);
    ec{k} = struct('reaction', tunedModel.ec.rxns{i}, 'kcat', num_value(tunedModel.ec.kcat(i)), ...
        'source', char(tunedModel.ec.source{i}), 'notes', char(tunedModel.ec.notes{i}));
end

out.raised = false;
out.tuned = tuned;
out.tuned_sorted = tuned(order);
out.final_growth = num_value(sol.f);
out.ec = ec;
end


function out = truncate_ties(inputs)
texts = inputs.truncate_inputs;
if ~iscell(texts)
    texts = cellstr(texts);
end
out = cell(1, numel(texts));
for i = 1:numel(texts)
    value = str2double(texts{i});
    try
        truncated = truncateValues({value}, 1);
        out{i} = struct('input', texts{i}, 'raised', false, 'output', num_value(truncated{1}));
    catch
        out{i} = struct('input', texts{i}, 'raised', true);
    end
end
end


function v = num_value(x)
% A double, or the text NaN, Inf or -Inf, which JSON cannot hold.
x = double(x);
if isnan(x)
    v = 'NaN';
elseif isinf(x)
    if x > 0
        v = 'Inf';
    else
        v = '-Inf';
    end
else
    v = x;
end
end


function restore_solver(previous)
if ~isempty(previous)
    setRavenSolver(previous);
end
end
