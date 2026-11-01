function results = flexibilize_edge_cases_ectestgem(ctx)
% MATLAB side of the flexibilize edge-case scenario; mirrors run.py. flexibilizeEnzConcs,
% relaxProteomicsGreedy and constrainEnzConcs on inputs the happy-path scenarios do not reach.
% A case that raises is a result, not a failure; flex_control returns real content.

previousSolver = getpref('RAVEN', 'solver', '');
restoreSolver = onCleanup(@() restore_solver(previousSolver));
setRavenSolver(char(ctx.inputs.matlab_solver));

inputs = ctx.inputs;
results = struct();

results.flex_control = run_case(@() flex_case(ctx, 'control'));
results.flex_growth_capped_by_uptake = run_case(@() flex_case(ctx, 'capped'));
results.flex_usage_bound_above_target = run_case(@() flex_case(ctx, 'preset'));
results.flex_infeasible_start = run_case(@() flex_case(ctx, 'infeasible'));
results.constrain_prior_lower_bound = run_case(@() prior_lb_case(ctx));
results.constrain_negative_concentration = run_case(@() negative_conc_case(ctx));

results.relax_control = run_case(@() relax_case(ctx, 'control'));
results.relax_tied_enzymes = run_case(@() relax_case(ctx, 'tied'));
results.relax_empty_enzyme_set = run_case(@() relax_case(ctx, 'empty'));

end


function out = flex_case(ctx, kind)
inputs = ctx.inputs;
expGrowth = double(inputs.exp_growth);
[m, ~, fixAdapter] = build_model(ctx, true);
switch kind
    case 'capped'
        % R1 carries the whole chain flux: cap growth at 5% above the enzyme-limited growth.
        s = solveLP(m);
        m.ub(strcmp(m.rxns, 'R1')) = 1.05 * s.f;
    case 'preset'
        conc = m.ec.concs(strcmp(m.ec.enzymes, 'P5'));
        m.ub(strcmp(m.rxns, 'usage_prot_P5')) = 10 * conc;
        expGrowth = double(inputs.exp_growth_high);
    case 'infeasible'
        m.lb(strcmp(m.rxns, 'R5')) = 50;
end
[mFlex, flexEnz] = flexibilizeEnzConcs(m, 'expGrowth', expGrowth, ...
    'foldChange', double(inputs.fold_change), 'iterPerEnzyme', double(inputs.iter_per_enzyme), ...
    'modelAdapter', fixAdapter, 'verbose', false);
n = numel(flexEnz.uniprotIDs);
records = cell(1, n);
for i = 1:n
    records{i} = struct('protein', flexEnz.uniprotIDs{i}, 'old_conc', num(flexEnz.oldConcs(i)), ...
        'flex_conc', num(flexEnz.flexConcs(i)), 'ratio_incr', num(flexEnz.ratioIncr(i)), ...
        'frequence', num(flexEnz.frequence(i)));
end
if n > 0
    [~, order] = sort(cellfun(@(r) r.protein, records, 'UniformOutput', false));
    records = records(order);
end
final = solveLP(mFlex);
optimal = final.stat == 1;
growth = 0;
if optimal
    growth = final.f;
end
out = struct('flexed', {records}, 'final_is_optimal', logical(optimal), ...
    'final_growth', num(growth), 'p5_usage_bounds', rxn_bounds(mFlex, 'usage_prot_P5'));
end


function out = prior_lb_case(ctx)
[m, ~] = build_model(ctx, false);
m.lb(strcmp(m.rxns, 'usage_prot_P5')) = 0.05;
m.ec.concs(strcmp(m.ec.enzymes, 'P5')) = 0.01;
m = constrainEnzConcs(m);
coeff = full(m.S(strcmp(m.mets, 'prot_pool'), strcmp(m.rxns, 'usage_prot_P5')));
out = struct('bounds', rxn_bounds(m, 'usage_prot_P5'), 'pool_coefficient', num(coeff));
end


function out = negative_conc_case(ctx)
[m, ~] = build_model(ctx, false);
m.ec.concs(strcmp(m.ec.enzymes, 'P4')) = -0.5;
m = constrainEnzConcs(m);
out = struct('bounds', rxn_bounds(m, 'usage_prot_P4'));
end


function out = relax_case(ctx, kind)
inputs = ctx.inputs;
[m, ~] = build_model(ctx, true);
target = double(inputs.relax_target);
args = {};
if strcmp(kind, 'tied')
    m = tie_p4(m);
    s = solveLP(m);
    target = 1.5 * s.f;
elseif strcmp(kind, 'empty')
    args = {'enzymeSet', {}};
end
result = relaxProteomicsGreedy(m, 'minimalGrowth', target, args{:});
relaxed = sort(fieldnames(result.relaxed));
nTrace = numel(result.trace);
if strcmp(kind, 'tied')
    % Which enzyme is relaxed first follows the dual solution of the LP, so only the
    % outcome of the loop is recorded.
    out = struct('relaxed', {reshape(relaxed, 1, [])}, 'final_growth', num(result.finalGrowth), ...
        'converged', logical(result.converged), 'n_steps', nTrace);
    return
end
trace = cell(1, nTrace);
for k = 1:nTrace
    t = result.trace(k);
    trace{k} = struct('iteration', num(t.iteration), 'relaxed', char(t.relaxedUniprot), ...
        'growth_before', num(t.growthBefore), 'growth_after', num(t.growthAfter), ...
        'shadow_price', num(t.shadowPrice));
end
out = struct('relaxed', {reshape(relaxed, 1, [])}, 'final_growth', num(result.finalGrowth), ...
    'converged', logical(result.converged), 'trace', {trace});
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
