model0 = pr_common('model'); res = struct();
revIdx = find(model0.rev == 1, 1); irrIdx = find(model0.rev == 0 & model0.lb == 0, 1);
mB = model0; mB.lb(revIdx) = 0; mB.ub(revIdx) = 1000;
mC = model0; mC.lb(irrIdx) = -1000;
res.variantRxns = struct('B_revWithLbZero', model0.rxns{revIdx}, 'C_irrevWithNegLb', model0.rxns{irrIdx});
vars = struct('A', model0, 'B', mB, 'C', mC);
for k = fieldnames(vars)'
    n = k{1}; m = vars.(n); r = struct();
    try
        evalc('[im, matchRev, r2i, i2r] = convertToIrrev(m);');
        r.rxnIds = im.rxns; r.matchRev = matchRev; r.nIrrev = numel(im.rxns);
        if n == 'B', sel = startsWith(im.rxns, model0.rxns{revIdx}); elseif n == 'C', sel = startsWith(im.rxns, model0.rxns{irrIdx}); else, sel = false(size(im.rxns)); end
        r.splitRxns = im.rxns(sel); r.splitLb = im.lb(sel); r.splitUb = im.ub(sel);
    catch e, r.error = pr_common('err', e); end
    res.(n) = r;
end
two = model0.rxns(find(model0.rev == 1, 2));
try
    evalc('[im, mr] = convertToIrrev(model0, two);');
    res.subset = struct('requested', {two}, 'rxnIds', {im.rxns}, 'matchRev', mr);
catch e, res.subset = struct('error', pr_common('err', e)); end
two2 = flipud(two);
try
    evalc('[im, mr] = convertToIrrev(model0, two2);');
    res.subsetReversedRequest = struct('requested', {two2}, 'rxnIds', {im.rxns}, 'matchRev', mr);
catch e, res.subsetReversedRequest = struct('error', pr_common('err', e)); end
pr_common('save', '2.6_convertToIrrev', res);
