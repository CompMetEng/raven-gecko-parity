function results = ec_lookup_edge_cases_ectestgem(ctx)
% MATLAB side of the EC-lookup and subsetting edge-case scenario; returns the shape run.py
% returns. A call that errors is reported as raised = true with an empty result; the message
% is not compared. Edits to the inputs sit outside the try blocks, so a mistake in them stops
% the scenario instead of looking like a refusal.

checks = ctx.inputs.database_checks;
for k = 1:numel(checks)
    results.database_checks.(checks(k).name) = database_check(ctx.inputs, checks(k));
end

results.apply_eccodes.full  = apply_eccodes(ctx.inputs, false);
results.apply_eccodes.light = apply_eccodes(ctx.inputs, true);

variants = cellstr(ctx.inputs.subset_variants);
for k = 1:numel(variants)
    results.subset.(variants{k}) = subset_check(ctx.inputs, variants{k});
end

end


function [adapter, model] = base(inputs)
adapter = ModelAdapterManager.getAdapter(inputs.adapter_matlab);
model = loadConventionalGEM('modelAdapter', adapter);
end


function out = database_check(inputs, check)
[adapter, model] = base(inputs);
ecModel = makeEcModel(model, false, adapter);
adapter.params.path = fullfile(inputs.variants_dir, check.variant);
if isempty(check.subset)
    mask = [];
else
    mask = ismember(ecModel.ec.rxns, cellstr(check.subset));
end
try
    ecModel = getECfromDatabase(ecModel, mask, check.action, adapter);
catch
    out.raised = true;
    out.eccodes = cell(1, 0);
    return
end
out.raised = false;
out.eccodes = eccode_records(ecModel);
end


function out = eccode_records(ecModel)
% ec.eccodes as one record per ec reaction, sorted by reaction id; an unassigned code is the
% empty string on both sides.
[rxnIds, order] = sort(ecModel.ec.rxns(:));
codes = ecModel.ec.eccodes(:);
codes = codes(order);
out = cell(1, numel(rxnIds));
for k = 1:numel(rxnIds)
    code = codes{k};
    if isempty(code)
        code = '';
    end
    out{k} = struct('reaction', rxnIds{k}, 'eccodes', code);
end
end


function out = apply_eccodes(inputs, light)
[adapter, model] = base(inputs);
ecModel = makeEcModel(model, light, adapter);
ecModel = getECfromGEM(ecModel);
try
    ecModel = applyECcodes(ecModel, false);
catch
    out.raised = true;
    out.reactions = cell(1, 0);
    return
end
[rxnIds, order] = sort(ecModel.rxns(:));
records = cell(1, numel(rxnIds));
for k = 1:numel(rxnIds)
    % Assigned by field: struct() with a cell value would build a struct array.
    record = struct('reaction', rxnIds{k});
    record.eccodes = reshape(eccode_tokens(ecModel, order(k)), 1, []);
    records{k} = record;
end
out.raised = false;
out.reactions = records;
end


function tokens = eccode_tokens(model, idx)
tokens = {};
if isfield(model, 'eccodes') && ~isempty(model.eccodes{idx})
    tokens = strsplit(model.eccodes{idx}, ';');
    tokens = sort(tokens(~cellfun(@isempty, tokens)));
end
end


function out = subset_check(inputs, variant)
[adapter, model] = base(inputs);
big = makeEcModel(model, false, adapter);
small = small_gem(variant, model);
try
    sub = getSubsetEcModel(big, small);
catch
    out.raised = true;
    out.model = struct();
    return
end
out.raised = false;
out.model = summary(sub);
end


function small = small_gem(variant, model)
switch variant
    case 'drop_reaction'
        small = removeReactions(model, {'R3'}, 'removeUnusedMets', true, 'removeUnusedGenes', true);
    case 'drop_gene'
        small = changeGrRules(model, {'R2'}, {'G1 and G2'}, true);
        nGenes = numel(small.genes);
        keep = ~strcmp(small.genes, 'G3');
        for f = {'genes', 'geneShortNames', 'geneMiriams', 'geneFrom', 'geneComps', 'proteins'}
            if isfield(small, f{1}) && size(small.(f{1}), 1) == nGenes
                small.(f{1}) = small.(f{1})(keep, :);
            end
        end
        small.rxnGeneMat = small.rxnGeneMat(:, keep);
    otherwise
        error('unknown subset variant %s', variant);
end
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
out.genes = reshape(sort(model.genes(:)), 1, []);
out.gene_associations = pairs(model.rxnGeneMat, model.rxns, model.genes, {'reaction', 'gene'}, false);
ec = model.ec;
out.ec_rxns    = reshape(sort(ec.rxns(:)), 1, []);
out.ec_enzymes = reshape(sort(ec.enzymes(:)), 1, []);
out.ec_genes   = reshape(sort(ec.genes(:)), 1, []);
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
