function results = enzyme_usage_edge_cases_ectestgem(ctx)
% MATLAB side of the enzyme-usage edge-case scenario; mirrors run.py. Every case builds its
% own ecModel (enzyme_usage_ectestgem's R2/R4-blocked single-route fixture) and records what
% the usage, report, bottleneck, pFBA and control-coefficient functions return for an input
% the happy-path scenario does not reach. A case that raises is a result, not a failure.

previousSolver = getpref('RAVEN', 'solver', '');
restoreSolver = onCleanup(@() restore_solver(previousSolver));
setRavenSolver(char(ctx.inputs.matlab_solver));

inputs = ctx.inputs;
results = struct();

[ecModel, ~] = build_model(ctx, true);
sol = solveLP(ecModel);
usage = getEnzymeUsage(ecModel, sol.x);
results.control_usage = run_case(@() struct('objective', num(sol.f), 'proteins', {usage_rows(usage)}));
results.control_report = run_case(@() report_struct(ecModel, usage, inputs, []));
results.control_bottlenecks = run_case(@() bottleneck_struct(ecModel, double(inputs.top_abs_usage)));
results.control_pfba = run_case(@() pfba_struct(ecModel, {}));

results.zero_ub_enzyme = run_case(@() zero_ub_case(ctx));
results.usage_exclude_zero = run_case(@() exclude_zero_case(ecModel, sol, inputs));
results.report_negative_top = run_case(@() report_struct(ecModel, usage, inputs, -1));
results.report_top_one = run_case(@() report_struct(ecModel, usage, inputs, 1));
results.usage_missing_reaction = run_case(@() missing_usage_case(ctx));

results.pfba_fraction_09 = run_case(@() pfba_struct(ecModel, {'fractionOfOptimum', 0.9}));
results.pfba_unknown_objective = run_case(@() pfba_struct(ecModel, {'rxnId', 'NOSUCH'}));
results.pfba_no_objective = run_case(@() pfba_no_objective_case(ctx));

results.tied_bottlenecks = run_case(@() tied_case(ctx));
results.control_coeffs_pool_binding = run_case(@() pool_binding_case(ctx));

end


function out = usage_rows(usageData)
n = numel(usageData.protID);
records = cell(1, n);
for k = 1:n
    records{k} = struct('protein', usageData.protID{k}, 'ub', num(usageData.UB(k)), ...
        'abs_usage', num(usageData.absUsage(k)), 'cap_usage', num(usageData.capUsage(k)));
end
if n > 0
    [~, order] = sort(usageData.protID);
    records = records(order);
end
out = records;
end


function out = table_rows(tbl, srcCol, dstField)
n = height(tbl);
records = cell(1, n);
for k = 1:n
    records{k} = struct('protein', tbl.protID{k}, 'abs_usage', num(tbl.absUsage(k)));
    records{k}.(dstField) = num(tbl.(srcCol)(k));
end
if n > 0
    [~, order] = sort(tbl.protID);
    records = records(order);
end
out = records;
end


function out = report_struct(ecModel, usageData, inputs, topOverride)
top = double(inputs.top_abs_usage);
if ~isempty(topOverride)
    top = topOverride;
end
report = reportEnzymeUsage(ecModel, usageData, ...
    'highCapUsage', double(inputs.high_cap_usage), 'topAbsUsage', top);
out = struct('high_cap_usage', {table_rows(report.highCapUsage, 'capUsage', 'cap_usage')}, ...
    'top_abs_usage', {table_rows(report.topAbsUsage, 'percUsage', 'perc_usage')}, ...
    'total_usage_flux', num(report.totalUsageFlux));
end


function out = bottleneck_struct(ecModel, top)
tbl = getEnzymeBottlenecks(ecModel, 'top', top);
n = height(tbl);
order = cell(1, n);
records = cell(1, n);
for k = 1:n
    order{k} = tbl.uniprot{k};
    records{k} = struct('protein', tbl.uniprot{k}, 'shadow_price', num(tbl.shadowPrice(k)), ...
        'flux', num(tbl.flux(k)), 'cap_usage', num(tbl.capUsage(k)), ...
        'upper_bound', num(tbl.upperBound(k)));
end
[~, sortIdx] = sort(order);
out = struct('order', {order}, 'rows', {records(sortIdx)});
end


function out = pfba_struct(ecModel, args)
optimum = solveLP(ecModel);
sol = getPfbaEnzymes(ecModel, args{:});
usageIdx = startsWith(ecModel.rxns, 'usage_prot_');
growth = sol.x(strcmp(ecModel.rxns, 'R5'));
out = struct('growth', num(growth), 'growth_over_optimum', num(growth / optimum.f), ...
    'usage_sum', num(sum(abs(sol.x(usageIdx)))));
end


function out = zero_ub_case(ctx)
[m, ~] = build_model(ctx, true);
m.ec.concs(strcmp(m.ec.enzymes, 'P1')) = 0;
m = constrainEnzConcs(m);
s = solveLP(m);
out = struct('proteins', {usage_rows(getEnzymeUsage(m, s.x))});
end


function out = exclude_zero_case(ecModel, sol, inputs)
u = getEnzymeUsage(ecModel, sol.x, 'zero', false);
out = struct('proteins', {usage_rows(u)}, 'report', report_struct(ecModel, u, inputs, []));
end


function out = missing_usage_case(ctx)
[m, ~] = build_model(ctx, true);
m = removeReactions(m, 'usage_prot_P1');
s = solveLP(m);
out = struct('proteins', {usage_rows(getEnzymeUsage(m, s.x))});
end


function out = pfba_no_objective_case(ctx)
[m, ~] = build_model(ctx, true);
m.c(:) = 0;
out = pfba_struct(m, {});
end


function out = tied_case(ctx)
[m, ~] = build_model(ctx, true);
m = tie_p4(m);
out = bottleneck_struct(m, double(ctx.inputs.top_abs_usage));
end


function out = pool_binding_case(ctx)
[m, ~] = build_model(ctx, false);
s = solveLP(m);
poolIdx = strcmp(m.rxns, 'prot_pool_exchange');
m.ub(poolIdx) = 0.5 * abs(s.x(poolIdx));
[enz, coeffs] = getConcControlCoeffs(m, 'foldChange', 2, 'limit', 0);
n = numel(m.ec.enzymes);
records = cell(1, n);
for k = 1:n
    records{k} = struct('protein', m.ec.enzymes{k}, 'analysed', logical(enz(k)), 'coeff', num(coeffs(k)));
end
[~, order] = sort(m.ec.enzymes);
out = struct('coeffs', {records(order)});
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
