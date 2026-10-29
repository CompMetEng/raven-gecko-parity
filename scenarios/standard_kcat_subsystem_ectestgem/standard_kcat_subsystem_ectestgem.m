function results = standard_kcat_subsystem_ectestgem(ctx)
% MATLAB side of the standard-kcat subsystem scenario; mirrors run.py. getStandardKcat on
% ecTestGEM plus a reaction without genes, with the subsystem assignments described in
% scenario.yml. A raise is a result.

adapter = ModelAdapterManager.getAdapter(ctx.inputs.adapter_matlab);
results.uniform = attempt(@() standard_case(adapter, struct('star', 'Alpha')));
results.case_variant = attempt(@() standard_case(adapter, struct('star', 'Alpha', 'R3', 'alpha')));
results.own_name = attempt(@() standard_case(adapter, struct('star', 'Alpha', 'R6', 'Beta')));
results.empty_subsystem = attempt(@() standard_case(adapter, struct('star', 'Alpha', 'R6', '')));
results.no_subsystem = attempt(@() standard_case(adapter, []));

end


function out = standard_case(adapter, subsystems)
model = loadConventionalGEM('modelAdapter', adapter);
rxnsToAdd.rxns      = {'R6'};
rxnsToAdd.equations = {'m1[c] => m2[c]'};
rxnsToAdd.grRules   = {''};
model = addRxns(model, rxnsToAdd, 3);
if ~isempty(subsystems)
    model.subSystems = cell(numel(model.rxns), 1);
    for i = 1:numel(model.rxns)
        name = subsystems.star;
        if isfield(subsystems, model.rxns{i})
            name = subsystems.(model.rxns{i});
        end
        if isempty(name)
            model.subSystems{i} = {};
        else
            model.subSystems{i} = {name};
        end
    end
end
ecModel = makeEcModel(model, false, adapter);
ecModel = getECfromGEM(ecModel);
ecModel.ec.kcat(:) = 10;
ecModel.ec.source(:) = {'manual'};
before = ecModel.ec.rxns(:);
ecModel = getStandardKcat(ecModel, 'modelAdapter', adapter);

n = numel(ecModel.ec.rxns);
rows = cell(1, n);
for i = 1:n
    src = '';
    if ~isempty(ecModel.ec.source{i})
        src = lower(ecModel.ec.source{i});
    end
    rows{i} = struct('reaction', ecModel.ec.rxns{i}, ...
        'kcat', struct('is_nan', logical(isnan(ecModel.ec.kcat(i))), 'value', nan2zero(ecModel.ec.kcat(i))), ...
        'source', src);
end
[~, order] = sort(ecModel.ec.rxns);
out = struct('new_rxns', {reshape(sort(setdiff(ecModel.ec.rxns(:), before)), 1, [])}, 'rows', {rows(order)});
end


function out = attempt(fn)
% A raise is a result on both sides.
try
    out = struct('raised', false, 'value', {fn()});
catch err
    fprintf('case raised: %s\n', err.message);
    out = struct('raised', true, 'value', struct());
end
end


function x = nan2zero(x)
x = double(x);
if isnan(x)
    x = 0;
end
end
