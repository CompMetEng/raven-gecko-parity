function results = ec_proteomics_ecyeastgem(ctx)
% MATLAB side of the genome-scale proteomics scenario; mirrors run.py. The stage-4 pipeline
% of GECKO's full_ecModel tutorial on ecYeastGEM, and loadProtData on the tutorial proteomics
% and on an edited file. Every case records the result of the call and a case that raises is
% a result. See scenario.yml.

previousSolver = getpref('RAVEN', 'solver', '');
restoreSolver = onCleanup(@() restore_solver(previousSolver));
setRavenSolver(char(ctx.inputs.matlab_solver));

inputs = ctx.inputs;
adapter = ModelAdapterManager.getAdapter(inputs.adapter_matlab);

results.loader = loader_cases(inputs, adapter);
results.enzyme_concs = enzyme_concs(inputs, adapter);
results.f_factor = f_factor_cases(inputs, adapter);
results.flux_data = flux_data_cases(inputs, adapter);
results.pipeline = run_case(@() pipeline_case(inputs, adapter));

end


function ecModel = fresh(inputs, adapter)
% The tutorial ecModel, growth as the objective.
ecModel = loadEcModel(char(inputs.model_file), adapter);
ecModel.c = double(strcmp(ecModel.rxns, char(inputs.bio_rxn)));
end


function out = solved(model)
% Whether the LP is feasible and its objective; an infeasible LP is a result.
sol = solveLP(model);
if isempty(sol.f) || sol.stat ~= 1
    out = struct('feasible', false, 'objective', num(0));
else
    out = struct('feasible', true, 'objective', num(sol.f));
end
end


function out = prot_rows(prot)
n = numel(prot.uniprotIDs);
nc = size(prot.abundances, 2);
rows = cell(1, n);
for i = 1:n
    row = cell(1, nc);
    for j = 1:nc
        row{j} = num(prot.abundances(i, j));
    end
    rows{i} = row;
end
out = struct('ids', {reshape(prot.uniprotIDs, 1, [])}, 'abundances', {rows});
end


function out = exchange_rows(model, inputs)
ids = cellstr(inputs.exchange_rxns);
out = cell(1, numel(ids));
for k = 1:numel(ids)
    idx = strcmp(model.rxns, ids{k});
    out{k} = struct('reaction', ids{k}, 'lower_bound', num(model.lb(idx)), 'upper_bound', num(model.ub(idx)));
end
end


function out = usage_bounds(model)
idx = find(startsWith(model.rxns, 'usage_prot_'));
[~, order] = sort(model.rxns(idx));
idx = idx(order);
out = cell(1, numel(idx));
for k = 1:numel(idx)
    out{k} = struct('reaction', model.rxns{idx(k)}, 'lower_bound', num(model.lb(idx(k))), ...
        'upper_bound', num(model.ub(idx(k))));
end
end


function out = concs_rows(model)
n = numel(model.ec.enzymes);
out = cell(1, n);
for k = 1:n
    out{k} = struct('protein', model.ec.enzymes{k}, 'conc', num(model.ec.concs(k)));
end
end


function out = loader_cases(inputs, adapter)
real = char(inputs.prot_data_file);
edge = char(inputs.prot_data_edge_file);
empty = char(inputs.prot_data_empty_file);
out.real_default = run_case(@() prot_rows(loadProtData(3, 'protDataFile', real, 'modelAdapter', adapter)));
out.real_unfiltered = run_case(@() prot_rows(loadProtData(3, 'protDataFile', real, 'filterData', false, 'modelAdapter', adapter)));
out.edge_default = run_case(@() prot_rows(loadProtData([3 2], 'protDataFile', edge, 'modelAdapter', adapter)));
out.edge_unfiltered = run_case(@() prot_rows(loadProtData([3 2], 'protDataFile', edge, 'filterData', false, 'modelAdapter', adapter)));
out.edge_max_missing_vector = run_case(@() prot_rows(loadProtData([3 2], 'protDataFile', edge, 'maxMissing', [0.5 1.0], 'modelAdapter', adapter)));
out.edge_strict = run_case(@() prot_rows(loadProtData([3 2], 'protDataFile', edge, 'minVal', 0.5, 'maxRSD', 0.3, 'cutLowest', 0, 'addStdevs', 0, 'modelAdapter', adapter)));
out.edge_cut_half = run_case(@() prot_rows(loadProtData([3 2], 'protDataFile', edge, 'cutLowest', 50, 'modelAdapter', adapter)));
out.edge_wrong_replicates = run_case(@() prot_rows(loadProtData([3 3], 'protDataFile', edge, 'modelAdapter', adapter)));
out.empty_cells_default = run_case(@() prot_rows(loadProtData([3 2], 'protDataFile', empty, 'modelAdapter', adapter)));
out.empty_cells_unfiltered = run_case(@() prot_rows(loadProtData([3 2], 'protDataFile', empty, 'filterData', false, 'modelAdapter', adapter)));
end


function out = enzyme_concs(inputs, adapter)
real = loadProtData(3, 'protDataFile', char(inputs.prot_data_file), 'modelAdapter', adapter);
out.real = run_case(@() concs_case(inputs, adapter, real, 1));
out.edge_condition_1 = run_case(@() concs_case(inputs, adapter, 'edge', 1));
out.edge_condition_2 = run_case(@() concs_case(inputs, adapter, 'edge', 2));
out.edge_out_of_range = run_case(@() concs_case(inputs, adapter, 'edge', 3));
end


function out = concs_case(inputs, adapter, prot, dataCol)
if ischar(prot)
    prot = loadProtData([3 2], 'protDataFile', char(inputs.prot_data_edge_file), 'modelAdapter', adapter);
end
model = fresh(inputs, adapter);
model = fillEnzConcs(model, prot, 'dataCol', dataCol);
concs = concs_rows(model);
model = constrainEnzConcs(model);
out = struct('n_measured', double(sum(~isnan(model.ec.concs))), 'concs', {concs}, ...
    'usage_bounds', {usage_bounds(model)}, 'growth', solved(model));
end


function out = f_factor_cases(inputs, adapter)
real = loadProtData(3, 'protDataFile', char(inputs.prot_data_file), 'modelAdapter', adapter);
model = fresh(inputs, adapter);
subset = model.ec.enzymes(1:100);
out.real = run_case(@() num(calculateFfactor(model, 'protData', real, 'modelAdapter', adapter)));
out.real_unfiltered = run_case(@() num(calculateFfactor(model, 'protData', ...
    loadProtData(3, 'protDataFile', char(inputs.prot_data_file), 'filterData', false, 'modelAdapter', adapter), ...
    'modelAdapter', adapter)));
out.edge = run_case(@() num(calculateFfactor(model, 'protData', ...
    loadProtData([3 2], 'protDataFile', char(inputs.prot_data_edge_file), 'modelAdapter', adapter), ...
    'modelAdapter', adapter)));
out.real_subset = run_case(@() num(calculateFfactor(model, 'protData', real, 'enzymes', subset, 'modelAdapter', adapter)));
end


function out = flux_data_cases(inputs, adapter)
data = loadFluxData('fluxDataFile', char(inputs.flux_data_file), 'modelAdapter', adapter);
n = size(data.exchFluxes, 1);
fluxes = cell(1, n);
for k = 1:n
    row = cell(1, size(data.exchFluxes, 2));
    for j = 1:numel(row)
        row{j} = num(data.exchFluxes(k, j));
    end
    fluxes{k} = row;
end
out.parsed = struct( ...
    'conds', {reshape(data.conds, 1, [])}, ...
    'p_tot', {num_list(data.Ptot)}, ...
    'gr_rate', {num_list(data.grRate)}, ...
    'exch_fluxes', {fluxes}, ...
    'exch_mets', {reshape(data.exchMets, 1, [])}, ...
    'exch_rxn_ids', {reshape(data.exchRxnIDs, 1, [])});
out.loose = run_case(@() constrained_case(inputs, adapter, data, 'loose'));
out.pct_10 = run_case(@() constrained_case(inputs, adapter, data, 10));
end


function out = num_list(values)
values = reshape(double(values), 1, []);
out = cell(1, numel(values));
for k = 1:numel(values)
    out{k} = num(values(k));
end
end


function out = constrained_case(inputs, adapter, data, looseStrict)
model = fresh(inputs, adapter);
model = constrainFluxData(model, 'fluxData', data, 'condition', 1, 'maxMinGrowth', 'max', ...
    'looseStrictFlux', looseStrict, 'modelAdapter', adapter);
out = struct('exchange', {exchange_rows(model, inputs)}, 'growth', solved(model));
end


function out = pipeline_case(inputs, adapter)
prot = loadProtData(3, 'protDataFile', char(inputs.prot_data_file), 'modelAdapter', adapter);
data = loadFluxData('fluxDataFile', char(inputs.flux_data_file), 'modelAdapter', adapter);
model = fresh(inputs, adapter);
model = fillEnzConcs(model, prot);
model = constrainEnzConcs(model);
f = calculateFfactor(model, 'protData', prot, 'modelAdapter', adapter);
model = setProtPoolSize(model, 'Ptot', data.Ptot(1), 'f', f, 'modelAdapter', adapter);
model = constrainFluxData(model, 'fluxData', data, 'condition', 1, 'maxMinGrowth', 'max', ...
    'looseStrictFlux', 'loose', 'modelAdapter', adapter);
poolIdx = strcmp(model.rxns, 'prot_pool_exchange');
out.f_factor = num(f);
out.pool_ub = num(model.ub(poolIdx));
out.growth_before = solved(model);
out.control_coeffs = run_case(@() ccc_rows(model, inputs));
out.flexibilize = run_case(@() flex_case(model, data, inputs, adapter));
end


function out = ccc_rows(model, inputs)
% The flag that says which enzymes were analysed follows the vertex of the LP; the
% coefficient is 0 for an enzyme that was not analysed.
[~, coeffs] = getConcControlCoeffs(model, 'foldChange', 2, 'limit', double(inputs.ccc_limit));
n = numel(model.ec.enzymes);
records = cell(1, n);
for k = 1:n
    records{k} = struct('protein', model.ec.enzymes{k}, 'coeff', num(coeffs(k)));
end
[~, order] = sort(model.ec.enzymes);
out = records(order);
end


function out = flex_case(model, data, inputs, adapter)
[flexModel, flexEnz] = flexibilizeEnzConcs(model, 'expGrowth', double(data.grRate(1)), ...
    'foldChange', double(inputs.fold_change), 'modelAdapter', adapter, 'verbose', false);
n = numel(flexEnz.uniprotIDs);
flexed = cell(1, n);
for i = 1:n
    flexed{i} = struct('protein', flexEnz.uniprotIDs{i}, ...
        'old_conc', num(flexEnz.oldConcs(i)), ...
        'flex_conc', num(flexEnz.flexConcs(i)), ...
        'ratio_incr', num(flexEnz.ratioIncr(i)), ...
        'frequence', double(flexEnz.frequence(i)));
end
poolIdx = strcmp(flexModel.rxns, 'prot_pool_exchange');
out = struct('flexed', {flexed}, 'usage_bounds', {usage_bounds(flexModel)}, ...
    'pool_ub', num(flexModel.ub(poolIdx)), 'growth_after', solved(flexModel));
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
