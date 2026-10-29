function results = ec_analysis_ecyeastgem(ctx)
% MATLAB side of the genome-scale ecModel analysis scenario; mirrors run.py. The post-solve
% analysis functions run on ecYeastGEM of GECKO's full_ecModel tutorial and every case
% records the result of the call; a case that raises is a result, not a failure. See
% scenario.yml for what is compared.

previousSolver = getpref('RAVEN', 'solver', '');
restoreSolver = onCleanup(@() restore_solver(previousSolver));
setRavenSolver(char(ctx.inputs.matlab_solver));

inputs = ctx.inputs;
adapter = ModelAdapterManager.getAdapter(ctx.inputs.adapter_matlab);
conv = loadConventionalGEM('modelAdapter', adapter);

results = struct();
ecModel = fresh(inputs, adapter);
sol = solveLP(ecModel);

% The usage functions and the mapping are given the same flux vector on both sides, read from
% data/growth_fluxes.tsv; the vector the LP returns differs between the two sides.
table = readtable(char(inputs.fluxes_file), 'FileType', 'text', 'Delimiter', '\t', 'TextType', 'char');
[found, loc] = ismember(ecModel.rxns, table.reaction);
shared = zeros(numel(ecModel.rxns), 1);
shared(found) = table.flux(loc(found));
usage = getEnzymeUsage(ecModel, shared);

results.growth = run_case(@() growth_struct(ecModel, sol, inputs));
results.usage = run_case(@() usage_rows(usage));
results.report = run_case(@() report_struct(ecModel, usage, inputs));
results.mapping = run_case(@() mapping_struct(ecModel, conv, shared));
solvedUsage = getEnzymeUsage(ecModel, sol.x);
results.usage_solved = run_case(@() usage_summary(solvedUsage));
results.bottlenecks = run_case(@() bottleneck_struct(ecModel, double(inputs.bottlenecks_top)));
results.pfba = run_case(@() pfba_struct(ecModel));
results.control_coeffs = run_case(@() control_struct(ecModel));
results.min_uptake = run_case(@() min_uptake_case(inputs, adapter));
results.fseof = run_case(@() fseof_case(inputs, adapter));
results.kcat_edits = run_case(@() kcat_edits_case(inputs, adapter));

end


function ecModel = fresh(inputs, adapter)
% The tutorial ecModel with glucose uptake open and growth as the objective.
ecModel = loadEcModel(char(inputs.model_file), adapter);
ecModel = set_bounds(ecModel, char(inputs.c_source), -1000, []);
ecModel = set_objective(ecModel, char(inputs.bio_rxn), 1);
end


function model = set_objective(model, rxnId, coefficient)
model.c = zeros(numel(model.rxns), 1);
model.c(strcmp(model.rxns, rxnId)) = coefficient;
end


function model = set_bounds(model, rxnId, lower, upper)
idx = strcmp(model.rxns, rxnId);
if ~isempty(lower) && lower > model.ub(idx)
    model.ub(idx) = lower;
end
if ~isempty(upper) && upper < model.lb(idx)
    model.lb(idx) = upper;
end
if ~isempty(lower)
    model.lb(idx) = lower;
end
if ~isempty(upper)
    model.ub(idx) = upper;
end
end


function out = growth_struct(ecModel, sol, inputs)
out = struct('objective', num(sol.f), 'exchange', {exchange_rows(ecModel, sol, inputs)});
end


function out = exchange_rows(ecModel, sol, inputs)
ids = cellstr(inputs.exchange_rxns);
out = cell(1, numel(ids));
for k = 1:numel(ids)
    out{k} = struct('reaction', ids{k}, 'flux', num(sol.x(strcmp(ecModel.rxns, ids{k}))));
end
end


function out = usage_summary(usageData)
% What every optimal vertex shares: the total usage and the largest usages as numbers.
values = sort(double(usageData.absUsage(:)), 'descend');
largest = cell(1, min(5, numel(values)));
for k = 1:numel(largest)
    largest{k} = num(values(k));
end
out = struct('total', num(sum(values)), 'largest', {largest});
end


function out = usage_rows(usageData)
n = numel(usageData.protID);
records = cell(1, n);
for k = 1:n
    records{k} = struct('protein', usageData.protID{k}, 'ub', num(usageData.UB(k)), ...
        'abs_usage', num(usageData.absUsage(k)), 'cap_usage', num(usageData.capUsage(k)));
end
[~, order] = sort(usageData.protID);
out = records(order);
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


function out = report_struct(ecModel, usageData, inputs)
report = reportEnzymeUsage(ecModel, usageData, ...
    'highCapUsage', double(inputs.high_cap_usage), 'topAbsUsage', double(inputs.top_abs_usage));
out = struct('high_cap_usage', {table_rows(report.highCapUsage, 'capUsage', 'cap_usage')}, ...
    'top_abs_usage', {table_rows(report.topAbsUsage, 'percUsage', 'perc_usage')}, ...
    'top_abs_usage_order', {reshape(report.topAbsUsage.protID, 1, [])}, ...
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


function out = pfba_struct(ecModel)
optimum = solveLP(ecModel);
sol = getPfbaEnzymes(ecModel);
usageIdx = startsWith(ecModel.rxns, 'usage_prot_');
objective = sum(ecModel.c .* sol.x);
out = struct('objective', num(objective), 'objective_over_optimum', num(objective / optimum.f), ...
    'usage_sum', num(sum(abs(sol.x(usageIdx)))));
end


function out = control_struct(ecModel)
[enz, coeffs] = getConcControlCoeffs(ecModel, 'foldChange', 2, 'limit', 0.001);
n = numel(ecModel.ec.enzymes);
records = cell(1, n);
for k = 1:n
    records{k} = struct('protein', ecModel.ec.enzymes{k}, 'analysed', logical(enz(k)), 'coeff', num(coeffs(k)));
end
[~, order] = sort(ecModel.ec.enzymes);
out = records(order);
end


function out = mapping_struct(ecModel, conv, fluxes)
[mappedFlux, enzUsageFlux, usageEnz] = mapRxnsToConv(ecModel, conv, fluxes);
n = numel(conv.rxns);
mapped = cell(1, n);
for k = 1:n
    mapped{k} = struct('reaction', conv.rxns{k}, 'flux', num(mappedFlux(k)));
end
m = numel(usageEnz);
usage = cell(1, m);
for k = 1:m
    usage{k} = struct('label', usageEnz{k}, 'flux', num(enzUsageFlux(k)));
end
[~, order] = sort(usageEnz);
out = struct('mapped', {mapped}, 'usage', {usage(order)});
end


function out = min_uptake_case(inputs, adapter)
m = fresh(inputs, adapter);
m = set_objective(m, char(inputs.c_source), 1);
m = set_bounds(m, char(inputs.bio_rxn), 0.25, []);
s = solveLP(m);
u = getEnzymeUsage(m, s.x);
out = struct('growth', growth_struct(m, s, inputs), 'usage', struct('total', usage_summary(u).total));
end


function out = fseof_case(inputs, adapter)
m = fresh(inputs, adapter);
fseof = ecFSEOF(m, char(inputs.fseof_target), char(inputs.c_source), ...
    'nSteps', double(inputs.fseof_steps), 'modelAdapter', adapter);
allTargets = [fseof.rxnTargets(:,1); fseof.transportTargets(:,1)];
levels = reshape(double(fseof.alpha), 1, []);
enforced = cell(1, numel(levels));
for k = 1:numel(levels)
    enforced{k} = num(levels(k));
end
out = struct('enforced_levels', {enforced}, ...
    'targets', {sort(reshape(allTargets, 1, []))}, ...
    'targets_ranked', {reshape(allTargets, 1, [])});
end


function picks = pick_edit_rows(ec)
% One rule, shared with run.py, from the order of ec.rxns.
rxns = ec.rxns;
perRow = full(sum(ec.rxnEnzMat ~= 0, 2));
kcat = double(ec.kcat(:));
iso = contains(rxns(:), '_EXP_');
rev = endsWith(rxns(:), '_REV');
picks.isozyme = first_of(iso & kcat > 0);
picks.reversible = first_of(rev & ~iso & kcat > 0);
picks.complex = first_of(perRow > 1 & kcat > 0 & ~iso & ~rev);
picks.zero = first_of(~iso & ~rev & perRow == 1 & kcat > 0);
end


function i = first_of(mask)
hits = find(mask);
if isempty(hits)
    i = 0;
else
    i = hits(1);
end
end


function out = column_struct(model, rxnId)
j = find(strcmp(model.rxns, rxnId));
[lowerValue, lowerKind] = bound(model.lb(j));
[upperValue, upperKind] = bound(model.ub(j));
rows = find(model.S(:, j));
ids = model.mets(rows);
[~, order] = sort(ids);
coefficients = cell(1, numel(rows));
for k = 1:numel(rows)
    r = rows(order(k));
    coefficients{k} = struct('metabolite', model.mets{r}, 'coefficient', double(model.S(r, j)));
end
out = struct('reaction', rxnId, 'lower_bound', lowerValue, 'lower_kind', lowerKind, ...
    'upper_bound', upperValue, 'upper_kind', upperKind, 'coefficients', {coefficients});
end


function [value, kind] = bound(x)
if isinf(x) && x > 0
    value = 0;
    kind = '+inf';
elseif isinf(x) && x < 0
    value = 0;
    kind = '-inf';
else
    value = double(x);
    kind = 'finite';
end
end


function out = kcat_edits_case(inputs, adapter)
m = fresh(inputs, adapter);
picks = pick_edit_rows(m.ec);
names = {'isozyme', 'reversible', 'complex', 'zero'};
factors = [2, 3, 0.5, 0];
picked = struct();
before = struct();
columnIds = {};
for k = 1:numel(names)
    i = picks.(names{k});
    if i > 0
        picked.(names{k}) = m.ec.rxns{i};
        before.(names{k}) = double(m.ec.kcat(i));
        m.ec.kcat(i) = m.ec.kcat(i) * factors(k);
        columnIds{end+1} = m.ec.rxns{i}; %#ok<AGROW>
    end
end
m = applyKcatConstraints(m);
sol = solveLP(m);
columns = cell(1, numel(columnIds));
for k = 1:numel(columnIds)
    columns{k} = column_struct(m, columnIds{k});
end
out = struct('picked', picked, 'kcat_before', before, 'columns', {columns}, 'growth', num(sol.f));
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
