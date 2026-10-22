model = pr_common('model'); G = model.genes(1:6);
ad.genes = G; ad.tissues = {'A';'B';'C'};
ad.levels = [0 0 0; NaN 5 5; 1 2 3; 2 2 2; 10 0 0; NaN NaN NaN]; ad.threshold = 1;
res = struct(); res.genes = G;
res.A = run1(model, ad, 'A', G);
res.B = run1(model, ad, 'B', G);
ad2 = ad; ad2.genes = [G; G(3)]; ad2.levels = [ad.levels; 3 2 1];
res.G3twice = run1(model, ad2, 'A', G);
ad3 = ad; ad3.levels(:,1) = NaN;
res.columnAAllNaN = run1(model, ad3, 'A', G);
res.hpa = struct('note', 'not run: no HPA fixture in tests/function_tests/test_data');
pr_common('save', '2.2_scoreModel', res);
function r = run1(model, ad, tissue, G)
r = struct();
try
    evalc('[rs, gs] = scoreModel(model, [], ad, tissue);');
    [~, ix] = ismember(G, model.genes);
    r.geneScoresG1toG6 = gs(ix); r.rxnScoresFirst20 = rs(1:min(20, numel(rs))); r.rxnScoresSum = sum(rs); r.nRxnScores = numel(rs);
    r.geneScoresAllUnique = unique(gs(:))';
catch e, r.error = pr_common('err', e); end
end
