function results = ec_tuning_ecyeastgem(ctx)
% MATLAB side of the genome-scale kcat tuning scenario; mirrors run.py. sensitivityTuning and
% sigmaFitter on the full_ecModel tutorial's ecYeastGEM before its tuning. See scenario.yml.

previousSolver = getpref('RAVEN', 'solver', '');
restoreSolver = onCleanup(@() restore_solver(previousSolver));
setRavenSolver(char(ctx.inputs.matlab_solver));

inputs = ctx.inputs;
adapter = ModelAdapterManager.getAdapter(inputs.adapter_matlab);

results = struct();
results.default = run_case(@() tune_case(inputs, adapter, double(inputs.desired_growth), [], {}));
if ~results.default.raised
    first = strsplit(results.default.value.tuned{1}.enzymes, ';');
    first = first{1};
    results.ignore_first = run_case(@() ignore_case(inputs, adapter, first));
else
    results.ignore_first = struct('raised', true, 'value', struct());
end
results.low_growth_fold_2 = run_case(@() tune_case(inputs, adapter, double(inputs.desired_growth_low), ...
    double(inputs.fold_change_low), {}));
results.sigma_fitter = run_case(@() sigma_case(inputs, adapter));

end


function ecModel = fresh(inputs, adapter)
% The model before tuning with glucose uptake open and the pool as the adapter sets it.
ecModel = loadEcModel(char(inputs.model_file), adapter);
ecModel.lb(strcmp(ecModel.rxns, char(inputs.c_source))) = -1000;
ecModel.c = double(strcmp(ecModel.rxns, char(inputs.bio_rxn)));
ecModel = setProtPoolSize(ecModel, [], [], [], adapter);
end


function out = tune_case(inputs, adapter, growth, fold, ignore)
model = fresh(inputs, adapter);
base = solveLP(model);
args = {'desiredGrowthRate', growth, 'modelAdapter', adapter, 'verbose', false};
if ~isempty(fold)
    args = [args, {'foldChange', fold}];
end
if ~isempty(ignore)
    args = [args, {'protToIgnore', ignore}];
end
[tunedModel, tunedKcats] = sensitivityTuning(model, args{:});
out = tuned_struct(tunedModel, tunedKcats);
out.base_growth = num(base.f);
end


function out = ignore_case(inputs, adapter, first)
out = tune_case(inputs, adapter, double(inputs.desired_growth), [], {first});
out.ignored = first;
end


function out = tuned_struct(model, tunedKcats)
n = numel(tunedKcats.rxns);
tuned = cell(1, n);
for i = 1:n
    tuned{i} = struct('reaction', tunedKcats.rxns{i}, 'enzymes', tunedKcats.enzymes{i}, ...
        'old_kcat', num(tunedKcats.oldKcat(i)), 'new_kcat', num(tunedKcats.newKcat(i)), ...
        'source', tunedKcats.source{i});
end
sol = solveLP(model);
out = struct('tuned', {tuned}, 'final_growth', num(sol.f), ...
    'kcat', {num_list(model.ec.kcat)}, ...
    'notes', {blanks_for_empty(model.ec.notes)}, ...
    'source', {blanks_for_empty(model.ec.source)});
end


function out = sigma_case(inputs, adapter)
model = fresh(inputs, adapter);
[fitted, sigma] = sigmaFitter(model, 'growthRate', double(inputs.sigma_growth_rate), ...
    'Ptot', double(inputs.sigma_p_tot), 'f', double(inputs.sigma_f), ...
    'makePlot', false, 'modelAdapter', adapter);
out = struct('sigma', num(sigma), ...
    'model_pool_ub', num(fitted.ub(strcmp(fitted.rxns, 'prot_pool_exchange'))));
end


function out = blanks_for_empty(values)
out = reshape(values, 1, []);
for k = 1:numel(out)
    if isempty(out{k})
        out{k} = '';
    end
end
end


function out = num_list(values)
values = reshape(double(values), 1, []);
out = cell(1, numel(values));
for k = 1:numel(values)
    out{k} = num(values(k));
end
end


function out = num(x)
% A number as a record, so that NaN is a value both sides can state.
x = double(x);
out = struct('is_nan', logical(isnan(x)), 'value', double(nan2zero(x)));
end


function x = nan2zero(x)
x(isnan(x)) = 0;
end


function out = run_case(fn)
% Result of one case; a raise is a result on both sides, not a failure.
try
    value = fn();
catch err
    fprintf('case raised: %s\n', err.message);
    out = struct('raised', true, 'value', struct());
    return
end
out = struct('raised', false, 'value', {value});
end


function restore_solver(previous)
if ~isempty(previous)
    setRavenSolver(previous);
end
end
