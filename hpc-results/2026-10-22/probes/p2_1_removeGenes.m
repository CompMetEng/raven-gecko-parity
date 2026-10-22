model = pr_common('model');
res = struct(); res.fieldsBefore = fieldnames(model); res.before = pr_common('stats', model);
v = {'default', {}; 'removeUnusedMets', {'removeUnusedMets', true}; 'removeBlockedRxns', {'removeBlockedRxns', true}; 'standardizeRulesFalse', {'standardizeRules', false}};
res.step1 = struct();
for i = 1:size(v,1)
    try
        evalc('m2 = removeGenes(model, {''YFR053C''}, v{i,2}{:});');
        res.step1.(v{i,1}) = pr_common('stats', m2);
    catch e, res.step1.(v{i,1}) = struct('error', pr_common('err', e)); end
end
% step 2: reactions with a single gene
nGenes = cellfun(@(r) numel(regexp(r, '[A-Za-z0-9_.]+(?= *(\)|$|and|or))', 'match')), model.grRules);
single = find(sum(model.rxnGeneMat ~= 0, 2) == 1);
res.step2 = struct('singleGeneRxnCount', numel(single));
if ~isempty(single)
    r = single(1); g = model.genes(model.rxnGeneMat(r,:) ~= 0);
    res.step2.rxn = model.rxns{r}; res.step2.gene = g; res.step2.rxnBefore = struct('lb', model.lb(r), 'ub', model.ub(r), 'grRule', model.grRules{r});
    try
        evalc('m2 = removeGenes(model, g);');
        k = find(strcmp(m2.rxns, model.rxns{r}));
        res.step2.rxnStays = ~isempty(k);
        if ~isempty(k), res.step2.rxnAfter = struct('lb', m2.lb(k), 'ub', m2.ub(k), 'grRule', m2.grRules{k}); end
        res.step2.after = pr_common('stats', m2);
    catch e, res.step2.error = pr_common('err', e); end
end
% step 3: gene used by no reaction
m3 = model; m3.genes{end+1,1} = 'UNUSED1';
m3.rxnGeneMat(:, end+1) = 0;
for f = fieldnames(m3)'
    fn = f{1};
    if startsWith(fn, 'gene') && ~strcmp(fn, 'genes') && numel(m3.(fn)) == numel(model.genes)
        if iscell(m3.(fn)), m3.(fn){end+1,1} = ''; else, m3.(fn)(end+1,1) = NaN; end
    end
end
res.step3 = struct();
try
    evalc('m4 = removeGenes(m3, {''UNUSED1''});');
    res.step3.unused1StillInGenes = any(strcmp(m4.genes, 'UNUSED1'));
    res.step3.genesAfter = numel(m4.genes);
    res.step3.fieldsAfter = fieldnames(m4);
    bad = {};
    for f = fieldnames(m4)'
        fn = f{1};
        if startsWith(fn, 'gene') && numel(m4.(fn)) ~= numel(m4.genes), bad{end+1,1} = fn; end %#ok<AGROW>
    end
    res.step3.geneFieldsWrongLength = bad;
catch e, res.step3.error = pr_common('err', e); end
res.step3.removeUsedGene.note = 'fieldnames before/after for the YFR053C removal';
try
    evalc('m5 = removeGenes(model, {''YFR053C''});');
    res.step4 = struct('fieldsBefore', {fieldnames(model)}, 'fieldsAfter', {fieldnames(m5)});
    bad = {};
    for f = fieldnames(m5)'
        fn = f{1};
        if startsWith(fn, 'gene') && numel(m5.(fn)) ~= numel(m5.genes), bad{end+1,1} = fn; end %#ok<AGROW>
    end
    res.step4.geneFieldsWrongLength = bad;
catch e, res.step4 = struct('error', pr_common('err', e)); end
pr_common('save', '2.1_removeGenes', res);
