function results = ec_model_build_light_humangem(ctx)
% MATLAB side of the light Human-GEM ecModel build scenario; mirrors run.py. The pipeline of
% GECKO's light_ecModel tutorial, from makeEcModel to setProtPoolSize, with a checkpoint
% after each stage. See scenario.yml.

previousSolver = getpref('RAVEN', 'solver', '');
restoreSolver = onCleanup(@() restore_solver(previousSolver));
setRavenSolver(char(ctx.inputs.matlab_solver));

adapter = ModelAdapterManager.getAdapter(ctx.inputs.adapter_matlab);
model = loadConventionalGEM('modelAdapter', adapter);

convCounts = struct('reactions', numel(model.rxns), 'metabolites', numel(model.mets), ...
    'genes', numel(model.genes));

ecModel = makeEcModel(model, true, adapter);
results.expansion = struct( ...
    'conv_counts', convCounts, ...
    'ec_counts', struct('ec_rxns', numel(ecModel.ec.rxns), 'enzymes', numel(ecModel.ec.enzymes), ...
        'genes', numel(ecModel.ec.genes)), ...
    'gecko_light', logical(ecModel.ec.geckoLight), ...
    'ec_rxns', {reshape(ecModel.ec.rxns, 1, [])}, ...
    'enzymes', {reshape(sort(ecModel.ec.enzymes(:)), 1, [])}, ...
    'n_model_reactions', numel(ecModel.rxns), ...
    'n_model_metabolites', numel(ecModel.mets));

ecModel = applyComplexData(ecModel, 'modelAdapter', adapter, 'verbose', false);
results.complexes = pairs(ecModel.ec.rxnEnzMat, ecModel.ec.rxns, ecModel.ec.enzymes, ...
    {'reaction', 'enzyme'}, true);

ecModel = getECfromGEM(ecModel);
codes = ecModel.ec.eccodes(:).';
for k = 1:numel(codes)
    if isempty(codes{k})
        codes{k} = '';
    end
end
results.eccodes = struct('populated', double(sum(~cellfun(@isempty, codes))), 'per_ec_rxn', {codes});

kcatList = fuzzyKcatMatching(ecModel, [], adapter);
dlkcat = readDLKcatOutput(ecModel, 'modelAdapter', adapter);
n = numel(dlkcat.rxns);
dlkcatRecords = cell(1, n);
for k = 1:n
    dlkcatRecords{k} = struct('reaction', dlkcat.rxns{k}, 'gene', dlkcat.genes{k}, 'kcat', num(dlkcat.kcats(k)));
end
results.dlkcat = sort_records(dlkcatRecords, {'reaction', 'gene'}, 'kcat');

merged = mergeDLKcatAndFuzzyKcats(dlkcat, kcatList);
m = numel(merged.rxns);
mergedRecords = cell(1, m);
for k = 1:m
    mergedRecords{k} = struct('reaction', merged.rxns{k}, 'kcat', num(merged.kcats(k)), ...
        'source', lower(merged.kcatSource{k}));
end

ecModel = assignKcatValues(ecModel, merged);
afterAssign = kcat_rows(ecModel);

% The standard-kcat step runs on a copy and its result is a checkpoint, whether it raises or
% not; the later stages continue from the model without it, so that they compare the same
% state on both sides.
standard = attempt(@() standard_step(ecModel, adapter));

ecModel = applyKcatConstraints(ecModel);
ecModel = setProtPoolSize(ecModel, [], [], [], adapter);

results.kcats = struct( ...
    'merged', {sort_records(mergedRecords, {'reaction', 'source'}, 'kcat')}, ...
    'after_assign', {afterAssign}, ...
    'standard', standard, ...
    'prot_coefficients', {prot_coefficients(ecModel)}, ...
    'nonzero_count', double(sum(ecModel.ec.kcat > 0)), ...
    'pool_ub', num(ecModel.ub(strcmp(ecModel.rxns, 'prot_pool_exchange'))));

end


function out = standard_step(ecModel, adapter)
before = ecModel.ec.rxns(:);
beforeEnzymes = ecModel.ec.enzymes(:);
probe = getStandardKcat(ecModel, 'modelAdapter', adapter);
out = struct('rows', {kcat_rows(probe)}, ...
    'new_rxns', {reshape(sort(setdiff(probe.ec.rxns(:), before)), 1, [])}, ...
    'new_enzymes', {reshape(sort(setdiff(probe.ec.enzymes(:), beforeEnzymes)), 1, [])});
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


function out = kcat_rows(ecModel)
n = numel(ecModel.ec.rxns);
rows = cell(1, n);
for i = 1:n
    src = '';
    if ~isempty(ecModel.ec.source{i})
        src = lower(ecModel.ec.source{i});
    end
    rows{i} = struct('reaction', ecModel.ec.rxns{i}, 'kcat', num(ecModel.ec.kcat(i)), 'source', src);
end
[~, order] = sort(ecModel.ec.rxns);
out = rows(order);
end


function out = prot_coefficients(model)
protIdx = find(startsWith(model.mets, 'prot_'));
records = {};
for p = 1:numel(protIdx)
    row = protIdx(p);
    [~, rxnIdx, vals] = find(model.S(row, :));
    for k = 1:numel(rxnIdx)
        records{end+1} = struct('reaction', model.rxns{rxnIdx(k)}, ...
            'metabolite', model.mets{row}, 'coefficient', full(double(vals(k)))); %#ok<AGROW>
    end
end
keys = cellfun(@(r) [r.reaction char(1) r.metabolite], records, 'UniformOutput', false);
[~, order] = sort(keys);
out = records(order);
end


function out = sort_records(records, keyFields, numField)
% Records sorted on the key fields joined with char(1), then on the number.
keys = cell(numel(records), 1);
for k = 1:numel(records)
    key = '';
    for f = 1:numel(keyFields)
        key = [key char(1) records{k}.(keyFields{f})]; %#ok<AGROW>
    end
    keys{k} = sprintf('%s%c%030.12f', key, char(1), records{k}.(numField).value);
end
[~, order] = sort(keys);
out = records(order);
end


function out = pairs(matrix, rowIds, colIds, names, withCoefficient)
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


function out = num(x)
% A number as a record, so that NaN is a value both sides can state.
x = double(x);
out = struct('is_nan', logical(isnan(x)), 'value', double(nan2zero(x)));
end


function x = nan2zero(x)
x(isnan(x)) = 0;
end


function restore_solver(previous)
if ~isempty(previous)
    setRavenSolver(previous);
end
end
