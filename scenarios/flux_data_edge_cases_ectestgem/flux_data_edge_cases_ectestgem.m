function results = flux_data_edge_cases_ectestgem(ctx)
% MATLAB side of the flux-data edge-case scenario; mirrors run.py. loadFluxData and
% constrainFluxData on files and arguments the happy-path scenario does not reach. A case that
% raises is a result, not a failure; load_control and constrain_control return real content.

previousSolver = getpref('RAVEN', 'solver', '');
restoreSolver = onCleanup(@() restore_solver(previousSolver));
setRavenSolver(char(ctx.inputs.matlab_solver));

results = struct();
results.load_control = run_case(@() load_case(ctx, 'control_two_cond.tsv'));
results.load_ragged_row = run_case(@() load_case(ctx, 'ragged.tsv'));
results.load_fixed_columns_only = run_case(@() load_case(ctx, 'fixed_only.tsv'));
results.load_nested_parentheses_header = run_case(@() load_case(ctx, 'paren_header.tsv'));
results.load_missing_file = run_case(@() load_case(ctx, 'does_not_exist.tsv'));

results.constrain_control = run_case(@() constrain_case(ctx, 'control_two_cond.tsv', 1, 'max'));
results.constrain_condition_by_name = run_case(@() constrain_case(ctx, 'control_two_cond.tsv', 'cond2', 'max'));
results.constrain_unknown_condition = run_case(@() constrain_case(ctx, 'control_two_cond.tsv', 'nope', 'max'));
results.constrain_missing_exchange = run_case(@() constrain_case(ctx, 'missing_exchange.tsv', 1, 'max'));
results.constrain_invalid_max_min_growth = run_case(@() constrain_case(ctx, 'control_two_cond.tsv', 1, 'median'));

end


function out = load_case(ctx, name)
adapter = ModelAdapterManager.getAdapter(ctx.inputs.adapter_matlab);
fluxData = loadFluxData('fluxDataFile', fullfile(char(ctx.inputs.data_dir), name), 'modelAdapter', adapter);
if isempty(fluxData)
    % loadFluxData warns and returns [] for a file it cannot read.
    out = struct('loaded', false);
    return
end
n = size(fluxData.exchFluxes, 1);
rows = cell(1, n);
for k = 1:n
    rows{k} = numlist(fluxData.exchFluxes(k, :));
end
out = struct('conds', {reshape(cellstr(fluxData.conds), 1, [])}, ...
    'p_tot', {numlist(fluxData.Ptot)}, 'gr_rate', {numlist(fluxData.grRate)}, ...
    'exch_fluxes', {rows}, 'exch_mets', {reshape(cellstr(fluxData.exchMets), 1, [])}, ...
    'exch_rxn_ids', {reshape(cellstr(fluxData.exchRxnIDs), 1, [])});
out.loaded = true;
end


function out = constrain_case(ctx, name, condition, maxMinGrowth)
adapter = ModelAdapterManager.getAdapter(ctx.inputs.adapter_matlab);
fluxData = loadFluxData('fluxDataFile', fullfile(char(ctx.inputs.data_dir), name), 'modelAdapter', adapter);
[m, ~, fixAdapter] = build_model(ctx, false);
m = constrainFluxData(m, 'fluxData', fluxData, 'condition', condition, ...
    'maxMinGrowth', maxMinGrowth, 'looseStrictFlux', 'loose', 'modelAdapter', fixAdapter);
out = struct('S1', rxn_bounds(m, 'S1'), 'S2', rxn_bounds(m, 'S2'), 'R5', rxn_bounds(m, 'R5'));
end


% --------------------------------------------------------------------------- %
% shared helpers (same in every probe scenario of this family)
% --------------------------------------------------------------------------- %

function out = num(x)
% A number as a record, so that NaN is a value both sides can state.
x = double(x);
out = struct('is_nan', logical(isnan(x)), 'value', double(nan2zero(x)));
end


function x = nan2zero(x)
x(isnan(x)) = 0;
end


function out = numlist(values)
values = reshape(double(values), 1, []);
out = cell(1, numel(values));
for k = 1:numel(values)
    out{k} = num(values(k));
end
end


function out = run_case(fn)
% Result of one case; a raise is a result on both sides, not a failure.
try
    value = fn();
catch
    out = struct('raised', true, 'value', struct());
    return
end
out = struct('raised', false, 'value', value);
end


function [ecModel, adapter, fixAdapter] = build_model(ctx, constrain, kcats)
% The R2/R4-blocked single-route ecTestGEM of enzyme_usage_ectestgem.
adapter = ModelAdapterManager.getAdapter(ctx.inputs.adapter_matlab);
fixAdapter = adapter;
fixAdapter.params.bioRxn = 'R5';
fixAdapter.params.c_source = 'S1';
model = loadConventionalGEM('modelAdapter', adapter);
for r = {'R2', 'R4'}
    idx = strcmp(model.rxns, r{1});
    model.lb(idx) = 0;
    model.ub(idx) = 0;
end
ecModel = makeEcModel(model, false, adapter);
ecModel = getECfromGEM(ecModel);
if nargin < 3 || isempty(kcats)
    kcats = ctx.inputs.kcats;
end
names = fieldnames(kcats);
for i = 1:numel(names)
    idx = strcmp(ecModel.ec.rxns, names{i});
    ecModel.ec.kcat(idx) = kcats.(names{i});
    ecModel.ec.source(idx) = {'manual'};
end
protData = loadProtData(1, [], [], adapter);
ecModel = fillEnzConcs(ecModel, protData);
if constrain
    ecModel = constrainEnzConcs(ecModel);
end
ecModel = applyKcatConstraints(ecModel);
ecModel = setProtPoolSize(ecModel, [], [], [], adapter);
end


function ecModel = tie_p4(ecModel)
% Make P4's cap equal to the flux it carries at the optimum, so that P4 and P5 both bind.
sol = solveLP(ecModel);
idx = strcmp(ecModel.rxns, 'usage_prot_P4');
ecModel.ec.concs(strcmp(ecModel.ec.enzymes, 'P4')) = abs(sol.x(idx));
ecModel = constrainEnzConcs(ecModel);
end


function out = rxn_bounds(model, rxnID)
idx = strcmp(model.rxns, rxnID);
out = struct('lb', num(model.lb(idx)), 'ub', num(model.ub(idx)));
end


function restore_solver(previous)
if ~isempty(previous)
    setRavenSolver(previous);
end
end
