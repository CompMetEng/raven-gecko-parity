model = pr_common('model'); res = struct();
% the nutrient inputs of smallYeast are closed (ub 0); open them so that the model grows
for q = {'glcIN','o2IN'}, model.ub(strcmp(model.rxns, q{1})) = 1000; end
res.openedInputs = {'glcIN','o2IN'};
res.biomassCandidates = model.rxns(model.c ~= 0);
bio = model.rxns{find(model.c ~= 0, 1)}; res.biomassRxn = bio;
sol = solveLP(model); res.growth = -sol.f;
neg = find(sol.x < -1e-6 & ~strcmp(model.rxns, bio)); res.negativeFluxRxns = model.rxns(neg);
if ~isempty(neg), target = model.rxns{neg(1)}; res.targetNote = 'carries negative flux at optimum';
else
    ex = getExchangeRxns(model); [~, ix] = ismember(ex, model.rxns); glc = model.rxns(ix(contains(lower(model.rxnNames(ix)), 'glucose')));
    res.glucoseExchangeCandidates = glc; target = glc{1}; res.targetNote = 'glucose exchange';
    mr = reverseRxns(model, target); res.reversedId = target; model = mr;
end
res.targetRxn = target;
res.step1 = fs(model, bio, target, 10, 0.9);
res.step2a = fs(model, bio, target, 10, 0.5);
res.step2b = fs(model, bio, target, 20, 0.9);
res.step2c = fs(model, bio, target, 20, 0.5);
mr = reverseRxns(model, target);
res.step3 = fs(mr, bio, target, 10, 0.9);
% Added to the requested steps: glcIN gives no target either way, so a product exchange that carries
% positive flux is used, once as it is and once reversed (flux then negative).
res.extra = struct();
base = model; if isfield(res, 'reversedId'), base = reverseRxns(model, res.reversedId); end
for q = {'ethOUT','acOUT','glyOUT'}
    t = q{1}; r = struct();
    sol = solveLP(base); r.growthWithoutTarget = -sol.f;
    r.asIs_ct09 = fs(base, bio, t, 10, 0.9); r.asIs_ct05 = fs(base, bio, t, 10, 0.5);
    br = reverseRxns(base, t);
    r.reversed_ct09 = fs(br, bio, t, 10, 0.9); r.reversed_ct05 = fs(br, bio, t, 10, 0.5);
    res.extra.(t) = r;
end
pr_common('save', '2.4_FSEOF', res);
function r = fs(model, bio, target, iters, ct)
r = struct('iterations', iters, 'corrThreshold', ct);
try
    evalc('t = FSEOF(model, bio, target, ''iterations'', iters, ''coefficient'', 0.9, ''corrThreshold'', ct);');
    r.fields = fieldnames(t); r.targetIds = model.rxns(t.logical); r.slope = t.slope(t.logical); r.knockdown = model.rxns(t.knockdown); r.knockout = model.rxns(t.knockout); r.amplify = model.rxns(t.amplify);
    r.allSlopesNonzero = sum(t.slope ~= 0);
    if isfield(t, 'corr'), r.corr = t.corr(t.logical); end
catch e, r.error = pr_common('err', e); end
end
