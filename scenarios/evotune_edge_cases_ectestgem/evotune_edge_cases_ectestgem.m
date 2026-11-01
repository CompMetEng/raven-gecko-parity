function results = evotune_edge_cases_ectestgem(ctx)
% MATLAB side of the evotune edge cases, mirroring run.py.

previousSolver = getpref('RAVEN', 'solver', '');
restoreSolver = onCleanup(@() restore_solver(previousSolver));
setRavenSolver(char(ctx.inputs.matlab_solver));

adapter = ModelAdapterManager.getAdapter(ctx.inputs.adapter_matlab);
model = loadConventionalGEM('modelAdapter', adapter);

inputs = ctx.inputs;
ev = inputs.evotune;

params = evotune_adapter(adapter, ev, ev.max_growth_weight);
heavy = evotune_adapter(adapter, ev, ev.heavy_max_growth_weight);
data = evo_data(inputs, 'flux_two_columns_file', {}, params);

results = struct();
[results.screen_control, ecModel, screen] = screen_case(model, adapter, inputs, inputs.control_sources, params, data);
results.screen_no_formulas = screen_case(model, adapter, inputs, inputs.control_sources, params, data, false);
results.screen_source_suffix = screen_case(model, adapter, inputs, inputs.suffix_sources, params, data);
results.screen_max_growth_weight = screen_case(model, adapter, inputs, inputs.control_sources, heavy, data);

zeroNames = {'known', 'unknown'};
zeroLists = {inputs.zero_flux_known, inputs.zero_flux_unknown};
for i = 1:numel(zeroNames)
    zeroData = evo_data(inputs, 'flux_one_column_file', zeroLists{i}, params);
    results.(['screen_zero_flux_' zeroNames{i}]) = screen_case(model, adapter, inputs, inputs.control_sources, params, zeroData);
end

shares = inputs.target_impact_shares;
for i = 1:numel(shares)
    key = ['select_mask_' strrep(sprintf('%.1f', shares(i)), '.', 'p')];
    mask = selectTunableMask(ecModel, screen, 'targetImpactShare', double(shares(i)));
    results.(key) = reshape(sort(ecModel.ec.rxns(mask)), 1, []);
end

flatModel = build_model(model, adapter, inputs, inputs.control_sources, true);
flatMask = ismember(flatModel.ec.rxns, cellstr(inputs.flat_mask));
try
    [~, tuned] = cmaesKcatTuning(flatModel, 'modelAdapter', params, 'evoData', data, ...
        'tunableMask', flatMask, 'seed', double(inputs.seed), 'verbose', false);
    flat.raised = false;
    flat.rxns = reshape(sort(tuned.rxns), 1, []);
    flat.n_generations = double(tuned.nGenerations);
    flat.converged = logical(tuned.converged);
    flat.kcat_unchanged = all(abs(tuned.newKcat - tuned.oldKcat) <= 1e-8 + 1e-5 * abs(tuned.oldKcat));
    results.cmaes_flat = flat;
catch
    results.cmaes_flat = struct('raised', true);
end

results.helpers = helpers(inputs);
end


function adapter = evotune_adapter(adapter, ev, maxGrowthWeight)
adapter.params.evotune = struct( ...
    'kcatSources', {reshape(cellstr(ev.source_groups), 1, [])}, ...
    'sigma0logSource', reshape(double(ev.sigma0_log_source), [], 1), ...
    'sigma0logDefault', double(ev.sigma0_log_default), ...
    'tieIsozymes', false, ...
    'maxGrowthWeight', double(maxGrowthWeight), ...
    'priorPenaltyWeight', 0, ...
    'maxGenerations', double(ev.max_generations), ...
    'rmseThreshold', double(ev.rmse_threshold));
end


function evoData = evo_data(inputs, fluxFileField, zeroFlux, adapter)
evoData.fluxData = loadFluxData(inputs.(fluxFileField), adapter);
evoData.fluxData.biomass = char(inputs.bio_rxn);
evoData.maxGrate = loadFluxData(inputs.max_growth_file, adapter);
evoData.maxGrate.biomass = char(inputs.bio_rxn);
evoData.zeroFlux = reshape(cellstr(zeroFlux), [], 1);
if isempty(zeroFlux)
    evoData.zeroFlux = {};
end
end


function ecModel = build_model(model, adapter, inputs, sources, withFormulas)
for id = {'R2', 'R4'}
    idx = strcmp(model.rxns, id{1});
    model.lb(idx) = 0; model.ub(idx) = 0;
end
ecModel = makeEcModel(model, false, adapter);
ecModel = getECfromGEM(ecModel);

kcatRxns = fieldnames(inputs.kcats);
for i = 1:numel(kcatRxns)
    idx = find(strcmp(ecModel.ec.rxns, kcatRxns{i}));
    if isempty(idx)
        error('no ec reaction %s', kcatRxns{i});
    end
    ecModel.ec.kcat(idx) = inputs.kcats.(kcatRxns{i});
end
sourceRxns = fieldnames(sources);
for i = 1:numel(sourceRxns)
    ecModel.ec.source{strcmp(ecModel.ec.rxns, sourceRxns{i})} = char(sources.(sourceRxns{i}));
end

if withFormulas
    % fillCarbonNum reads model.metFormulas, which ecTestGEM does not have; empty formulas give every exchange a carbon weight of 1.
    ecModel.metFormulas = repmat({''}, numel(ecModel.mets), 1);
end
ecModel.c = double(strcmp(ecModel.rxns, char(inputs.bio_rxn)));
ecModel = applyKcatConstraints(ecModel);
pool = inputs.pool;
ecModel = setProtPoolSize(ecModel, double(pool.p_tot), double(pool.f), double(pool.sigma), adapter);
end


function [out, ecModel, screen] = screen_case(model, adapter, inputs, sources, params, data, withFormulas)
if nargin < 7
    withFormulas = true;
end
ecModel = build_model(model, adapter, inputs, sources, withFormulas);
screen = [];
try
    screen = screenKcatLeverage(ecModel, 'modelAdapter', params, 'evoData', data);
catch
    out = struct('raised', true);
    return
end

n = height(screen);
rows = cell(1, n);
for i = 1:n
    rows{i} = struct('rxn', screen.rxnId{i}, 'n_isozymes', double(screen.nIsozymes(i)), ...
        'source_group', screen.sourceGroup{i}, 'kcat0', screen.kcat0(i), ...
        'leverage', screen.leverage(i), 'sigma0_log', screen.sigma0log(i), ...
        'rank', screen.rank(i), 'cum_leverage_share', screen.cumLeverageShare(i));
end
[~, order] = sort(cellfun(@(r) r.rxn, rows, 'UniformOutput', false));
out.raised = false;
out.order = reshape(screen.rxnId, 1, []);
out.rows = rows(order);
end


function out = helpers(inputs)
kcat = reshape(double(inputs.helper_kcat), [], 1);
kcat0 = reshape(double(inputs.helper_kcat0), [], 1);
sigma = reshape(double(inputs.helper_sigma0_log), [], 1);

scoreFcn = @(vec) mean(abs(log(vec(:) ./ kcat)));
points = parsimonyFrontier(scoreFcn, kcat, kcat0, sigma);
best = bestParsimonious(points);

tieRxns = cellstr(string(inputs.helper_tie_rxns));
tieKcat0 = reshape(double(inputs.helper_tie_kcat0), [], 1);
tieSources = cellstr(string(inputs.helper_tie_sources));
tieMap.without_sources = reshape(double(isozymeTieMap(tieRxns, tieKcat0)) - 1, 1, []);
tieMap.with_sources = reshape(double(isozymeTieMap(tieRxns, tieKcat0, 'sources', tieSources)) - 1, 1, []);
try
    isozymeTieMap(tieRxns, tieKcat0(1:end-1));
    tieMap.length_mismatch = 'not raised';
catch
    tieMap.length_mismatch = 'raised';
end
out.tie_map = tieMap;

out.n_changed = double(nChanged(kcat, kcat0));
frontier = cell(1, numel(points));
for i = 1:numel(points)
    frontier{i} = struct('threshold', double(points(i).threshold), 'n_changed', double(points(i).nChanged), ...
        'mean_dev', points(i).meanDev, 'rmse', points(i).rmse, 'kcat', reshape(points(i).kcat, 1, []));
end
out.frontier = frontier;
out.best_threshold = double(best.threshold);

rmse = reshape(double(inputs.helper_points_rmse), [], 1);
nCh = reshape(double(inputs.helper_points_n_changed), [], 1);
rmse(double(inputs.helper_nan_index)) = NaN;
nanPoints = repmat(struct('threshold', [], 'kcat', [], 'nChanged', [], 'meanDev', [], 'rmse', []), numel(rmse), 1);
for i = 1:numel(rmse)
    nanPoints(i) = struct('threshold', i - 1, 'kcat', kcat0, 'nChanged', nCh(i), 'meanDev', 0, 'rmse', rmse(i));
end
try
    chosen = bestParsimonious(nanPoints);
    out.best_with_nan = struct('raised', false, 'threshold', double(chosen.threshold));
catch
    out.best_with_nan = struct('raised', true);
end
end


function restore_solver(previous)
if ~isempty(previous)
    setRavenSolver(previous);
end
end
