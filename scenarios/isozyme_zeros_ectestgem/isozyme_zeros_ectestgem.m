function results = isozyme_zeros_ectestgem(ctx)
% MATLAB side of the measured-zero and isozyme-coverage scenario; mirrors run.py.
% getIsozymeCoverage and zeroUnmeasuredIsozymes with ec.concs written per case (NaN not
% measured, 0 measured absent), constrainEnzConcs and flexibilizeEnzConcs with a measured zero,
% and levelComplexConcs with a zero subunit. A case that raises is a result, not a failure.

previousSolver = getpref('RAVEN', 'solver', '');
restoreSolver = onCleanup(@() restore_solver(previousSolver));
setRavenSolver(char(ctx.inputs.matlab_solver));

inputs = ctx.inputs;
limit = double(inputs.detection_limit);
results = struct();

results.coverage_partial = run_case(@() coverage_record(build_model(ctx, 'partial', true, '')));
results.coverage_measured_zero_entry = run_case(@() coverage_record(build_model(ctx, 'measured_zero_entry', true, '')));
results.coverage_zero_subunit = run_case(@() coverage_record(build_model(ctx, 'zero_subunit', true, '')));
results.coverage_all_measured = run_case(@() coverage_record(build_model(ctx, 'all_measured', true, '')));
results.coverage_none_measured = run_case(@() coverage_record(build_model(ctx, 'none_measured', true, '')));
results.coverage_shared = run_case(@() coverage_record(build_model(ctx, 'partial', true, 'shared')));

results.zero_method_zero = run_case(@() zero_record(build_model(ctx, 'partial', true, ''), 'method', 'zero'));
results.zero_detection_limit_default = run_case(@() zero_record(build_model(ctx, 'partial', true, ''), 'method', 'detectionLimit'));
results.zero_detection_limit_explicit = run_case(@() zero_record(build_model(ctx, 'partial', true, ''), 'method', 'detectionLimit', 'limit', limit));
results.zero_bound_not_raised = run_case(@() not_raised(ctx, limit));
results.zero_shared_enzyme_unchanged = run_case(@() zero_record(build_model(ctx, 'partial', true, 'shared'), 'method', 'zero'));
results.zero_no_partial_group = run_case(@() zero_record(build_model(ctx, 'all_measured', true, ''), 'method', 'zero'));
results.zero_unknown_method = run_case(@() zero_record(build_model(ctx, 'partial', true, ''), 'method', 'bogus'));
results.zero_negative_limit = run_case(@() zero_record(build_model(ctx, 'partial', true, ''), 'method', 'detectionLimit', 'limit', -1));
results.zero_no_concentration_no_limit = run_case(@() zero_record(build_model(ctx, 'none_measured', true, ''), 'method', 'detectionLimit'));

results.constrain_measured_zero = run_case(@() constrain_zero(ctx));
results.flexibilize_measured_zero = run_case(@() flex_zero(ctx));
results.level_zero_subunit = run_case(@() level_zero(ctx));

end


function out = coverage_record(m)
[coverage, zeroable] = getIsozymeCoverage(m);
n = height(coverage);
rows = cell(1, n);
for i = 1:n
    rows{i} = struct('rxn', char(coverage.rxn{i}), 'n_entries', double(coverage.nEntries(i)), ...
        'n_measured', double(coverage.nMeasured(i)), 'n_unmeasured', double(coverage.nUnmeasured(i)), ...
        'measured', {sorted_row(coverage.measuredEnzymes{i})}, ...
        'unmeasured', {sorted_row(coverage.unmeasuredEnzymes{i})}, ...
        'shared', {sorted_row(coverage.sharedEnzymes{i})});
end
if n > 0
    [~, order] = sort(cellfun(@(r) r.rxn, rows, 'UniformOutput', false));
    rows = rows(order);
end
out = struct('rows', {rows}, 'zeroable', {sorted_row(zeroable)});
end


function out = zero_record(m, varargin)
[m, changed] = zeroUnmeasuredIsozymes(m, varargin{:});
n = height(changed);
records = cell(1, n);
for i = 1:n
    records{i} = struct('uniprot', char(changed.uniprot{i}), 'new_ub', num(changed.newUpperBound(i)));
end
if n > 0
    [~, order] = sort(cellfun(@(r) r.uniprot, records, 'UniformOutput', false));
    records = records(order);
end
out = struct('changed', {records}, 'usage_ub', {usage_ub_list(m)}, 'concs', {concs_list(m)});
end


function out = not_raised(ctx, limit)
m = build_model(ctx, 'partial', true, '');
m.ub(strcmp(m.rxns, 'usage_prot_P3')) = 0.05;
out = zero_record(m, 'method', 'detectionLimit', 'limit', limit);
end


function out = constrain_zero(ctx)
m = build_model(ctx, 'zero_p5', true, '');
out = struct('usage_ub', {usage_ub_list(m)}, 'growth', growth(m));
end


function out = flex_zero(ctx)
inputs = ctx.inputs;
[m, ~, fixAdapter] = build_model(ctx, 'zero_p5', true, '');
lastwarn('');
[m, flexEnz] = flexibilizeEnzConcs(m, 'expGrowth', double(inputs.exp_growth), ...
    'foldChange', double(inputs.fold_change), 'iterPerEnzyme', double(inputs.iter_per_enzyme), ...
    'modelAdapter', fixAdapter, 'verbose', false);
[msg, id] = lastwarn;
listed = {};
if strcmp(id, 'flexibilizeEnzConcs:measuredZero')
    tok = regexp(msg, 'enzymes (.*?) have a concentration of 0', 'tokens', 'once', 'ignorecase');
    if ~isempty(tok)
        listed = strsplit(strrep(tok{1}, ' ', ''), ',');
    end
end
out = struct('warned_enzymes', {sorted_row(listed)}, 'flexed', {sorted_row(flexEnz.uniprotIDs)}, ...
    'p5_ub', num(m.ub(strcmp(m.rxns, 'usage_prot_P5'))), 'final_growth', growth(m));
end


function out = level_zero(ctx)
m = build_model(ctx, 'zero_subunit', false, 'pair');
[m, result] = levelComplexConcs(m);
out = struct('concs', {concs_list(m)}, 'adjusted', {sorted_row(result.adjusted)}, ...
    'added', {sorted_row(result.added)});
end


function g = growth(m)
sol = solveLP(m);
value = 0;
if sol.stat == 1
    value = sol.f;
end
g = num(abs(value));
end


function out = sorted_row(c)
if isempty(c)
    out = {};
    return
end
out = reshape(sort(cellstr(c)), 1, []);
end


function out = concs_list(m)
proteins = {'P1', 'P2', 'P3', 'P4', 'P5'};
out = cell(1, numel(proteins));
for i = 1:numel(proteins)
    out{i} = struct('protein', proteins{i}, 'conc', num(m.ec.concs(strcmp(m.ec.enzymes, proteins{i}))));
end
end


function out = usage_ub_list(m)
proteins = {'P1', 'P2', 'P3', 'P4', 'P5'};
out = cell(1, numel(proteins));
for i = 1:numel(proteins)
    out{i} = struct('protein', proteins{i}, 'ub', num(m.ub(strcmp(m.rxns, ['usage_prot_' proteins{i}]))));
end
end


function m = set_layout(m, layout)
% Replace the rows of ec.rxnEnzMat by the copies of the layout.
matrix = zeros(size(m.ec.rxnEnzMat));
rxnNames = fieldnames(layout);
for i = 1:numel(rxnNames)
    row = find(strcmp(m.ec.rxns, rxnNames{i}));
    subunits = layout.(rxnNames{i});
    enzNames = fieldnames(subunits);
    for j = 1:numel(enzNames)
        matrix(row, strcmp(m.ec.enzymes, enzNames{j})) = double(subunits.(enzNames{j}));
    end
end
m.ec.rxnEnzMat = sparse(matrix);
end


% --------------------------------------------------------------------------- %
% shared helpers (same in every probe scenario of this family)
% --------------------------------------------------------------------------- %

function out = num(x)
% A number as a record, so that NaN is a value both sides can state.
x = double(x);
out = struct('is_nan', logical(isnan(x)), 'value', double(nan2zero(x)));
end


function x = nan2zero(x)
x(isnan(x)) = 0;
end


function out = run_case(fn)
% Result of one case; a raise is a result on both sides, not a failure.
try
    value = fn();
catch
    out = struct('raised', true, 'value', struct());
    return
end
out = struct('raised', false, 'value', value);
end


function [ecModel, adapter, fixAdapter] = build_model(ctx, concsName, constrain, layoutName)
% The R2/R4-blocked single-route ecTestGEM with the concentrations of the case.
adapter = ModelAdapterManager.getAdapter(ctx.inputs.adapter_matlab);
fixAdapter = adapter;
fixAdapter.params.bioRxn = 'R5';
fixAdapter.params.c_source = 'S1';
model = loadConventionalGEM('modelAdapter', adapter);
for r = {'R2', 'R4'}
    idx = strcmp(model.rxns, r{1});
    model.lb(idx) = 0;
    model.ub(idx) = 0;
end
ecModel = makeEcModel(model, false, adapter);
ecModel = getECfromGEM(ecModel);
kcats = ctx.inputs.kcats;
names = fieldnames(kcats);
for i = 1:numel(names)
    idx = strcmp(ecModel.ec.rxns, names{i});
    ecModel.ec.kcat(idx) = kcats.(names{i});
    ecModel.ec.source(idx) = {'manual'};
end
protData = loadProtData(1, [], [], adapter);
ecModel = fillEnzConcs(ecModel, protData);
if ~isempty(layoutName)
    ecModel = set_layout(ecModel, ctx.inputs.layouts.(layoutName));
end
concs = ctx.inputs.concs.(concsName);
proteins = fieldnames(concs);
for i = 1:numel(proteins)
    value = concs.(proteins{i});
    if isempty(value)
        value = NaN;
    end
    ecModel.ec.concs(strcmp(ecModel.ec.enzymes, proteins{i})) = double(value);
end
if constrain
    ecModel = constrainEnzConcs(ecModel);
end
ecModel = applyKcatConstraints(ecModel);
ecModel = setProtPoolSize(ecModel, [], [], [], adapter);
end


function restore_solver(previous)
if ~isempty(previous)
    setRavenSolver(previous);
end
end
