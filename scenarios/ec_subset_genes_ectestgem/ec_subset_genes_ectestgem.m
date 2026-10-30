function results = ec_subset_genes_ectestgem(ctx)
% MATLAB side of the ecModel subset scenario; mirrors run.py. getSubsetEcModel on ecTestGEM
% with context models that lack a reaction and a gene of its rule. A raise is a result. See
% scenario.yml.

adapter = ModelAdapterManager.getAdapter(ctx.inputs.adapter_matlab);
model = loadConventionalGEM('modelAdapter', adapter);

variants = {'identity', 'drop_reaction', 'drop_reaction_and_gene'};
flavours = {'full', 'light'};
for f = 1:numel(flavours)
    big = makeEcModel(model, strcmp(flavours{f}, 'light'), adapter);
    for i = 1:numel(variants)
        small = small_gem(variants{i}, model);
        try
            sub = getSubsetEcModel(big, small);
            results.(flavours{f}).(variants{i}) = struct('raised', false, 'model', summary(sub));
        catch err
            fprintf('case raised: %s\n', err.message);
            results.(flavours{f}).(variants{i}) = struct('raised', true, 'model', struct());
        end
    end
end

end


function small = small_gem(variant, model)
switch variant
    case 'identity'
        small = model;
    case 'drop_reaction'
        small = removeReactions(model, {'R2'}, 'removeUnusedMets', true, 'removeUnusedGenes', true);
    case 'drop_reaction_and_gene'
        small = removeReactions(model, {'R2'}, 'removeUnusedMets', true, 'removeUnusedGenes', false);
        nGenes = numel(small.genes);
        keep = ~strcmp(small.genes, 'G3');
        for f = {'genes', 'geneShortNames', 'geneMiriams', 'geneFrom', 'geneComps', 'proteins'}
            if isfield(small, f{1}) && size(small.(f{1}), 1) == nGenes
                small.(f{1}) = small.(f{1})(keep, :);
            end
        end
        small.rxnGeneMat = small.rxnGeneMat(:, keep);
end
end


function out = summary(model)
assoc = {};
[rows, cols] = find(model.rxnGeneMat);
for k = 1:numel(rows)
    assoc{end+1} = [model.rxns{rows(k)} ':' model.genes{cols(k)}]; %#ok<AGROW>
end
out = struct( ...
    'reactions', {reshape(sort(model.rxns), 1, [])}, ...
    'genes', {reshape(sort(model.genes), 1, [])}, ...
    'gene_associations', {reshape(sort(assoc), 1, [])}, ...
    'ec_rxns', {reshape(sort(model.ec.rxns), 1, [])}, ...
    'ec_genes', {reshape(sort(model.ec.genes), 1, [])}, ...
    'ec_enzymes', {reshape(sort(model.ec.enzymes), 1, [])});
end
