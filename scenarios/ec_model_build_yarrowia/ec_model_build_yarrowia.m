function results = ec_model_build_yarrowia(ctx)
% MATLAB side of the Yarrowia ecModel build scenario; mirrors run.py. The ecModel pipeline on
% iYali5-GEM with the data of the eciYali5-GEM repository, with a checkpoint after each stage.
% See scenario.yml.

previousSolver = getpref('RAVEN', 'solver', '');
restoreSolver = onCleanup(@() restore_solver(previousSolver));
setRavenSolver(char(ctx.inputs.matlab_solver));

[adapterFolder, adapterClass] = fileparts(char(ctx.inputs.adapter_matlab));
if ~contains(path, adapterFolder)
    addpath(adapterFolder);
end
root = char(ctx.inputs.eciyali_dir);

results.as_is = pipeline(ctx, adapterClass, root);

% The UniProt and KEGG snapshots quote the cells that hold several values, such as
% "1.1.1.25; 2.5.1.19" in the EC column. `unquoted` runs the pipeline on a copy of the data
% folder with every double quote removed, which separates what depends on the quotes from what
% does not.
copy = tempname;
mkdir(fullfile(copy, 'data'));
cleanup = onCleanup(@() rmdir(copy, 's'));
for name = {'uniprot.tsv', 'kegg.tsv', 'customKcats.tsv'}
    text = fileread(fullfile(root, 'data', name{1}));
    fid = fopen(fullfile(copy, 'data', name{1}), 'w', 'n', 'UTF-8');
    fprintf(fid, '%s', strrep(text, '"', ''));
    fclose(fid);
end
results.unquoted = pipeline(ctx, adapterClass, copy);

end


function results = pipeline(ctx, adapterClass, dataRoot)
% The pipeline with data/ read from dataRoot.
adapter = feval(adapterClass, char(ctx.inputs.eciyali_dir));
adapter.params.path = dataRoot;
model = loadConventionalGEM('modelAdapter', adapter);

convCounts = struct('reactions', numel(model.rxns), 'metabolites', numel(model.mets), ...
    'genes', numel(model.genes));

ecModel = makeEcModel(model, false, adapter);
results.expansion = struct( ...
    'conv_counts', convCounts, ...
    'ec_counts', struct('ec_rxns', numel(ecModel.ec.rxns), 'enzymes', numel(ecModel.ec.enzymes), ...
        'genes', numel(ecModel.ec.genes)), ...
    'ec_rxns', {reshape(ecModel.ec.rxns, 1, [])}, ...
    'enzymes', {reshape(ecModel.ec.enzymes, 1, [])}, ...
    'mw', {num_list(ecModel.ec.mw)}, ...
    'n_model_reactions', numel(ecModel.rxns), ...
    'n_model_metabolites', numel(ecModel.mets));

% The model has no EC codes of its own. getECfromGEM raises for that; whether the call raises
% is a result, and the codes after it are recorded either way.
fromGemRaised = false;
try
    ecModel = getECfromGEM(ecModel);
catch err
    fprintf('case raised: %s\n', err.message);
    fromGemRaised = true;
end
fromGem = eccode_list(ecModel);
ecModel = getECfromDatabase(ecModel, [], 'display', adapter);
fromDatabase = eccode_list(ecModel);
results.eccodes = struct('from_gem_raised', fromGemRaised, 'from_gem', {fromGem}, 'from_database', {fromDatabase}, ...
    'populated_after_database', double(sum(~cellfun(@isempty, fromDatabase))));

kcatList = fuzzyKcatMatching(ecModel, [], adapter);
n = numel(kcatList.rxns);
fuzzyRecords = cell(1, n);
for k = 1:n
    code = '';
    if ~isempty(kcatList.eccodes{k})
        code = char(kcatList.eccodes{k});
    end
    origin = -1;
    if ~isnan(kcatList.origin(k))
        origin = double(kcatList.origin(k));
    end
    wildcard = -1;
    if ~isnan(kcatList.wildcardLvl(k))
        wildcard = double(kcatList.wildcardLvl(k));
    end
    fuzzyRecords{k} = struct('reaction', kcatList.rxns{k}, 'kcat', num(kcatList.kcats(k)), ...
        'eccode', code, 'origin', origin, 'wildcard_level', wildcard);
end
[~, fuzzyOrder] = sort(kcatList.rxns);
fuzzyRecords = fuzzyRecords(fuzzyOrder);

ecModel = assignKcatValues(ecModel, kcatList);
afterAssign = kcat_rows(ecModel);

ecModel = applyCustomKcats(ecModel, 'modelAdapter', adapter);
afterCustom = kcat_rows(ecModel);

% The standard-kcat step runs on a copy and its result is a checkpoint, whether it raises or
% not; the later stages continue from the model without it, so that they compare the same
% state on both sides.
standard = attempt(@() standard_step(ecModel, adapter));

ecModel = applyKcatConstraints(ecModel);
ecModel = setProtPoolSize(ecModel, [], [], [], adapter);

results.kcats = struct( ...
    'fuzzy', {fuzzyRecords}, ...
    'after_assign', {afterAssign}, ...
    'after_custom', {afterCustom}, ...
    'standard', standard, ...
    'prot_coefficients', {prot_coefficients(ecModel)}, ...
    'nonzero_count', double(sum(ecModel.ec.kcat > 0)), ...
    'pool_ub', num(ecModel.ub(strcmp(ecModel.rxns, 'prot_pool_exchange'))));

end


function out = eccode_list(ecModel)
out = ecModel.ec.eccodes(:).';
for k = 1:numel(out)
    if isempty(out{k})
        out{k} = '';
    end
end
end


function out = num_list(values)
values = reshape(double(values), 1, []);
out = cell(1, numel(values));
for k = 1:numel(values)
    out{k} = num(values(k));
end
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
        src = ecModel.ec.source{i};
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
