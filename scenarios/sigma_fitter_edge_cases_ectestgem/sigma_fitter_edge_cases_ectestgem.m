function results = sigma_fitter_edge_cases_ectestgem(ctx)
% MATLAB side of the sigma-fit edge cases, mirroring run.py. sigmaFitter has one algorithm, so each
% case is run once and returned under both the _default and the _grid key.

previousSolver = getpref('RAVEN', 'solver', '');
restoreSolver = onCleanup(@() restore_solver(previousSolver));
setRavenSolver(char(ctx.inputs.matlab_solver));

adapter = ModelAdapterManager.getAdapter(ctx.inputs.adapter_matlab);
model = loadConventionalGEM('modelAdapter', adapter);

inputs = ctx.inputs;
uptake = struct('S1', [double(inputs.uptake_lower_bound), 1000]);
fixed = struct(char(inputs.bio_rxn), [double(inputs.fixed_biomass), double(inputs.fixed_biomass)]);

names = {'on_grid', 'between_grid', 'uptake_cap', 'fixed_biomass'};
bounds = {struct(), struct(), uptake, fixed};
growth = {inputs.growth_on_grid, inputs.growth_between_grid_points, ...
    inputs.growth_above_uptake_limit, inputs.growth_between_grid_points};

results = struct();
for i = 1:numel(names)
    out = run_case(model, adapter, inputs, bounds{i}, double(growth{i}), strcmp(names{i}, 'fixed_biomass'));
    results.([names{i} '_default']) = out;
    results.([names{i} '_grid']) = out;
end
end


function out = run_case(model, adapter, inputs, bounds, growthRate, reportBounds)
% Single-route ecModel (R2 and R4 blocked) with optional bounds edited on the conventional GEM first.
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

ecModel = makeEcModel(model, false, adapter);
ecModel = getECfromGEM(ecModel);
kcatRxns = fieldnames(inputs.kcats);
for i = 1:numel(kcatRxns)
    idx = strcmp(ecModel.ec.rxns, kcatRxns{i});
    ecModel.ec.kcat(idx) = inputs.kcats.(kcatRxns{i});
    ecModel.ec.source(idx) = {'manual'};
end
ecModel.c = double(strcmp(ecModel.rxns, char(inputs.bio_rxn)));
ecModel = applyKcatConstraints(ecModel);
ecModel = setProtPoolSize(ecModel, [], [], [], adapter);

[fittedModel, sigma] = sigmaFitter(ecModel, ...
    'growthRate', growthRate, ...
    'Ptot', double(inputs.p_tot), 'f', double(inputs.f), ...
    'makePlot', false, 'modelAdapter', adapter);

poolIdx = strcmp(fittedModel.rxns, 'prot_pool_exchange');
bioIdx = strcmp(fittedModel.rxns, char(inputs.bio_rxn));
out.sigma = double(sigma);
out.model_pool_ub = double(fittedModel.ub(poolIdx));
out.objective_fixed = fittedModel.lb(bioIdx) == fittedModel.ub(bioIdx);
if reportBounds
    out.objective_lb = double(fittedModel.lb(bioIdx));
    out.objective_ub = double(fittedModel.ub(bioIdx));
end
end


function restore_solver(previous)
if ~isempty(previous)
    setRavenSolver(previous);
end
end
