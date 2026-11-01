function results = ec_fva_ecyeastgem(ctx)
% MATLAB side of the genome-scale ecFVA scenario; mirrors run.py. ecFVA on ecYeastGEM
% constrained by the tutorial's flux data, once with growth free below the measured rate and
% once with growth fixed at it. See scenario.yml.
%
% ecFVA starts a parallel pool of the job's cores if none is running.

previousSolver = getpref('RAVEN', 'solver', '');
restoreSolver = onCleanup(@() restore_solver(previousSolver));
setRavenSolver(char(ctx.inputs.matlab_solver));

inputs = ctx.inputs;
adapter = ModelAdapterManager.getAdapter(inputs.adapter_matlab);
conv = loadConventionalGEM('modelAdapter', adapter);

fluxData = loadFluxData('fluxDataFile', char(inputs.flux_data_file), 'modelAdapter', adapter);
fluxData.grRate(1) = double(inputs.growth_rate);

ecModel = loadEcModel(char(inputs.model_file), adapter);
ecModel = constrainFluxData(ecModel, 'fluxData', fluxData, 'condition', 1, ...
    'maxMinGrowth', 'max', 'looseStrictFlux', 'loose', 'modelAdapter', adapter);
bioIdx = strcmp(ecModel.rxns, char(inputs.bio_rxn));

results.free = fva_case(ecModel, conv, bioIdx);

ecModel.lb(bioIdx) = ecModel.ub(bioIdx);
results.fixed = fva_case(ecModel, conv, bioIdx);

end


function out = fva_case(ecModel, conv, bioIdx)
[minFlux, maxFlux] = ecFVA(ecModel, conv);
n = numel(conv.rxns);
rows = cell(1, n);
for k = 1:n
    rows{k} = struct('reaction', conv.rxns{k}, 'min_flux', num(minFlux(k)), 'max_flux', num(maxFlux(k)));
end
out = struct('growth_bounds', {{num(ecModel.lb(bioIdx)), num(ecModel.ub(bioIdx))}}, 'fva', {rows});
end


function out = num(x)
% A number as a record, so that NaN is a value both sides can state.
x = double(x);
if isnan(x)
    x = 0;
    out = struct('is_nan', true, 'value', x);
else
    out = struct('is_nan', false, 'value', x);
end
end


function restore_solver(previous)
if ~isempty(previous)
    setRavenSolver(previous);
end
end
