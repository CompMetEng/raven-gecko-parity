function results = ec_model_edit_edge_cases_ectestgem(ctx)
% MATLAB side of the ecModel edit edge-case scenario; returns the shape run.py returns.
% A call that errors is reported as raised = true with an empty result; the message is not
% compared. The edits to the inputs sit outside the try blocks, so a mistake in them stops
% the scenario instead of looking like a refusal.

% The solver is named by the scenario; the preference is global and the nightly runs several
% scenarios in one MATLAB session, so whatever was set is put back on the way out.
previousSolver = getpref('RAVEN', 'solver', '');
restoreSolver = onCleanup(@() restore_solver(previousSolver));
setRavenSolver(char(ctx.inputs.matlab_solver));

cases = ctx.inputs.add_cases;
for k = 1:numel(cases)
    results.add_new_rxns.(cases(k).name) = add_case(ctx.inputs, cases(k));
end

idCases = ctx.inputs.id_cases;
for k = 1:numel(idCases)
    adapter = ModelAdapterManager.getAdapter(ctx.inputs.adapter_matlab);
    model = loadConventionalGEM('modelAdapter', adapter);
    model.rxns{strcmp(model.rxns, idCases(k).from)} = idCases(k).to;
    results.ids_with_suffix_text.(idCases(k).name).map_rxns = map_rxns(model, adapter);
    results.ids_with_suffix_text.(idCases(k).name).ec_fva   = ec_fva_case(ctx.inputs, model, adapter);
end

end


function restore_solver(previous)
if ~isempty(previous)
    setRavenSolver(previous);
end
end


function out = add_case(inputs, c)
adapter = ModelAdapterManager.getAdapter(inputs.adapter_matlab);
model = loadConventionalGEM('modelAdapter', adapter);
ecModel = makeEcModel(model, false, adapter);

spec = inputs.new_reaction;
newRxns.rxns      = {spec.id};
newRxns.rxnNames  = {spec.name};
newRxns.equations = {spec.equation};
newRxns.grRules   = {c.gr_rule};

% jsondecode turns an array of objects into a struct array and an empty array into [].
enzSpecs = c.enzymes;
if isempty(enzSpecs)
    newEnzymes.enzymes = cell(0, 1);
    newEnzymes.genes   = cell(0, 1);
    newEnzymes.mw      = zeros(0, 1);
else
    newEnzymes.enzymes = arrayfun(@(e) e.enzyme, enzSpecs, 'UniformOutput', false);
    newEnzymes.genes   = arrayfun(@(e) e.gene,   enzSpecs, 'UniformOutput', false);
    newEnzymes.mw      = arrayfun(@(e) e.mw,     enzSpecs);
    newEnzymes.enzymes = newEnzymes.enzymes(:);
    newEnzymes.genes   = newEnzymes.genes(:);
    newEnzymes.mw      = newEnzymes.mw(:);
end

try
    [ecModel, rxnsAdded, enzAdded] = addNewRxnsToEC(ecModel, newRxns, newEnzymes, adapter);
catch
    out.raised = true;
    out.rxns_added = cell(1, 0);
    out.enz_added = cell(1, 0);
    out.model = struct();
    return
end
out.raised = false;
out.rxns_added = reshape(sort(rxnsAdded(:)), 1, []);
out.enz_added  = reshape(sort(enzAdded(:)), 1, []);
out.model = summary(ecModel);
end


function out = map_rxns(model, adapter)
ecModel = makeEcModel(model, false, adapter);
[sortedIds, ~] = sort(ecModel.rxns(:));
rank = containers.Map(sortedIds, num2cell(1:numel(sortedIds)));
fluxVect = zeros(numel(ecModel.rxns), 1);
for k = 1:numel(ecModel.rxns)
    fluxVect(k) = rank(ecModel.rxns{k});
end
try
    [mappedFlux, enzUsageFlux, usageEnz] = mapRxnsToConv(ecModel, model, fluxVect);
catch
    out.raised = true;
    out.mapped = cell(1, 0);
    out.usage = cell(1, 0);
    return
end
[convIds, order] = sort(model.rxns(:));
mappedFlux = mappedFlux(order);
mapped = cell(1, numel(convIds));
for k = 1:numel(convIds)
    mapped{k} = struct('reaction', convIds{k}, 'flux', double(mappedFlux(k)));
end
[usageIds, order] = sort(usageEnz(:));
enzUsageFlux = enzUsageFlux(order);
usage = cell(1, numel(usageIds));
for k = 1:numel(usageIds)
    usage{k} = struct('enzyme', usageIds{k}, 'flux', double(enzUsageFlux(k)));
end
out.raised = false;
out.mapped = mapped;
out.usage = usage;
end


function out = ec_fva_case(inputs, model, adapter)
ecModel = makeEcModel(model, false, adapter);
ecModel = getECfromGEM(ecModel);
kcats = inputs.kcats;
kcatRxns = fieldnames(kcats);
for i = 1:numel(kcatRxns)
    idx = strcmp(ecModel.ec.rxns, kcatRxns{i});
    ecModel.ec.kcat(idx) = kcats.(kcatRxns{i});
    ecModel.ec.source(idx) = {'manual'};
end
ecModel = applyKcatConstraints(ecModel);
ecModel = setProtPoolSize(ecModel, [], [], [], adapter);
try
    [minFlux, maxFlux] = ecFVA(ecModel, model);
catch
    out.raised = true;
    out.fva = cell(1, 0);
    return
end
n = numel(model.rxns);
records = cell(1, n);
for k = 1:n
    records{k} = struct('reaction', model.rxns{k}, ...
        'min_flux', double(minFlux(k)), 'max_flux', double(maxFlux(k)));
end
[~, order] = sort(model.rxns);
out.raised = false;
out.fva = records(order);
end


function out = summary(model)
[rxnIds, order] = sort(model.rxns(:));
reactions = cell(1, numel(rxnIds));
for k = 1:numel(rxnIds)
    i = order(k);
    [lowerValue, lowerKind] = bound(model.lb(i));
    [upperValue, upperKind] = bound(model.ub(i));
    reactions{k} = struct( ...
        'id',          rxnIds{k}, ...
        'lower_bound', lowerValue, ...
        'lower_kind',  lowerKind, ...
        'upper_bound', upperValue, ...
        'upper_kind',  upperKind);
end
out.reactions = reactions;
out.metabolites = reshape(sort(model.mets(:)), 1, []);
out.gene_associations = pairs(model.rxnGeneMat, model.rxns, model.genes, {'reaction', 'gene'}, false);
ec = model.ec;
out.ec_rxns = reshape(sort(ec.rxns(:)), 1, []);
[enzIds, order] = sort(ec.enzymes(:));
genes = ec.genes(:);
mw = double(ec.mw(:));
enzymes = cell(1, numel(enzIds));
for k = 1:numel(enzIds)
    i = order(k);
    enzymes{k} = struct('enzyme', enzIds{k}, 'gene', genes{i}, 'mw', mw(i));
end
out.ec_enzymes = enzymes;
out.ec_coupling = pairs(ec.rxnEnzMat, ec.rxns, ec.enzymes, {'reaction', 'enzyme'}, true);
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


function out = pairs(matrix, rowIds, colIds, names, withCoefficient)
% Non-zero entries as records sorted on the two keys joined with char(1), which reproduces
% Python's tuple ordering.
[rows, cols, values] = find(matrix);
keys = cell(numel(rows), 1);
for k = 1:numel(rows)
    keys{k} = [rowIds{rows(k)} char(1) colIds{cols(k)}];
end
[~, order] = sort(keys);
out = cell(1, numel(order));
for k = 1:numel(order)
    i = order(k);
    record = struct(names{1}, rowIds{rows(i)}, names{2}, colIds{cols(i)});
    if withCoefficient
        record.coefficient = double(values(i));
    end
    out{k} = record;
end
end
