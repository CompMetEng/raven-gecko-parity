function results = complex_levelling_ectestgem(ctx)
% MATLAB side of the complex-levelling scenario; mirrors run.py. getComplexFamilies,
% levelComplexConcs and flexibilizeEnzConcs with perComplex on ecTestGEM with complexes written
% into ec.rxnEnzMat. A case that raises is a result, not a failure.

previousSolver = getpref('RAVEN', 'solver', '');
restoreSolver = onCleanup(@() restore_solver(previousSolver));
setRavenSolver(char(ctx.inputs.matlab_solver));

inputs = ctx.inputs;
results = struct();

layouts = {'pair', 'trio', 'shared', 'isozyme', 'binding'};
for k = 1:numel(layouts)
    layout = layouts{k};
    results.(['families_' layout]) = run_case(@() families_case(ctx, layout));
end

for layout = {'pair', 'trio'}
    for statistic = {'median', 'mean'}
        name = ['level_' layout{1} '_' statistic{1}];
        results.(name) = run_case(@() level_case(ctx, layout{1}, {}, 'statistic', statistic{1}));
    end
end

results.level_shared_subunit = run_case(@() level_case(ctx, 'shared', {}));
results.level_isozyme_variant = run_case(@() level_case(ctx, 'isozyme', {}));
results.level_unmeasured_fill = run_case(@() level_case(ctx, 'trio', {'P3'}, 'fillMissing', true));
results.level_unmeasured_nofill = run_case(@() level_case(ctx, 'trio', {'P3'}, 'fillMissing', false));
results.level_complex_unmeasured = run_case(@() level_case(ctx, 'trio', {'P1', 'P2', 'P3'}));
results.level_unknown_statistic = run_case(@() level_case(ctx, 'trio', {}, 'statistic', 'mode'));

targets = {'measured', 'exp_growth'; 'median', 'exp_growth'; 'mean', 'exp_growth_mean'};
for t = 1:size(targets, 1)
    label = targets{t, 1};
    target = double(inputs.(targets{t, 2}));
    if strcmp(label, 'measured')
        statistic = '';
    else
        statistic = label;
    end
    for perComplex = [false true]
        if perComplex
            suffix = 'per_complex';
        else
            suffix = 'default';
        end
        name = ['flex_' label '_' suffix];
        results.(name) = run_case(@() flex_case(ctx, statistic, target, perComplex));
    end
end

end


function out = families_case(ctx, layout)
[m, ~] = build_model(ctx);
m = set_layout(m, ctx.inputs.layouts.(layout));
out = families_record(getComplexFamilies(m), m);
end


function out = level_case(ctx, layout, nanProteins, varargin)
[m, ~] = build_model(ctx);
m = set_layout(m, ctx.inputs.layouts.(layout));
for i = 1:numel(nanProteins)
    m.ec.concs(strcmp(m.ec.enzymes, nanProteins{i})) = NaN;
end
before = concs_list(m);
[m, result] = levelComplexConcs(m, varargin{:});
after = concs_list(m);
m = finish_model(ctx, m);
skipped = cellfun(@(c) strjoin(sort(c(:))', char(1)), result.skipped, 'UniformOutput', false);
skipped = sort(skipped);
skippedLists = cell(1, numel(skipped));
for i = 1:numel(skipped)
    skippedLists{i} = reshape(strsplit(skipped{i}, char(1)), 1, []);
end
out = struct('concs_before', {before}, 'concs', {after}, 'usage_ub', {usage_ub_list(m)}, ...
    'adjusted', {sorted_row(result.adjusted)}, 'added', {sorted_row(result.added)}, ...
    'shared', {sorted_row(result.shared)}, 'skipped', {skippedLists});
end


function out = flex_case(ctx, statistic, expGrowth, perComplex)
% An empty statistic leaves the measured concentrations of the complex as they are.
inputs = ctx.inputs;
[m, ~, fixAdapter] = build_model(ctx);
m = set_layout(m, inputs.layouts.binding);
if ~isempty(statistic)
    m = levelComplexConcs(m, 'statistic', statistic);
end
m = finish_model(ctx, m);
start = solveLP(m);
[mFlex, flexEnz] = flexibilizeEnzConcs(m, 'expGrowth', expGrowth, ...
    'foldChange', double(inputs.fold_change), 'iterPerEnzyme', double(inputs.iter_per_enzyme), ...
    'perComplex', perComplex, 'modelAdapter', fixAdapter, 'verbose', false);
ids = reshape(flexEnz.uniprotIDs, 1, []);
n = numel(ids);
records = cell(1, n);
for i = 1:n
    records{i} = struct('protein', ids{i}, 'old_conc', num(flexEnz.oldConcs(i)), ...
        'flex_conc', num(flexEnz.flexConcs(i)), 'ratio_incr', num(flexEnz.ratioIncr(i)), ...
        'frequence', num(flexEnz.frequence(i)));
end
if n > 0
    [~, order] = sort(cellfun(@(r) r.protein, records, 'UniformOutput', false));
    records = records(order);
end
final = solveLP(mFlex);
optimal = final.stat == 1;
growth = 0;
if optimal
    growth = final.f;
end
out = struct('concs', {concs_list(m)}, 'start_growth', num(start.f), 'flexed', {records}, ...
    'flex_set', {sorted_row(ids)}, 'flex_order', {ids}, 'final_is_optimal', logical(optimal), ...
    'final_growth', num(growth), 'usage_ub', {usage_ub_list(mFlex)});
end


function out = families_record(families, m)
nC = numel(families.complexes);
complexes = cell(1, nC);
keys = cell(1, nC);
for c = 1:nC
    cplx = families.complexes{c};
    [sortedEnz, order] = sort(cplx.enzymes(:));
    subunits = cell(1, numel(order));
    for i = 1:numel(order)
        subunits{i} = struct('enzyme', sortedEnz{i}, 'copies', num(cplx.copies(order(i))), ...
            'owned', logical(cplx.owned(order(i))));
    end
    complexes{c} = struct('subunits', {subunits});
    keys{c} = strjoin(sortedEnz', char(1));
end
if nC > 0
    [~, order] = sort(keys);
    complexes = complexes(order);
end
enzymes = sort(m.ec.enzymes(:))';
members = cell(1, numel(enzymes));
for i = 1:numel(enzymes)
    members{i} = struct('enzyme', enzymes{i}, 'members', '');
    [~, e] = ismember(enzymes{i}, m.ec.enzymes);
    c = families.owner(e);
    if c > 0
        cplx = families.complexes{c};
        members{i}.members = strjoin(sort(cplx.enzymes(cplx.owned))', char(1));
    end
end
out = struct('complexes', {complexes}, 'owner_members', {members}, ...
    'shared', {sorted_row(families.shared)});
end


function out = sorted_row(values)
% A sorted 1xN cell array of char, also for an empty input.
if isempty(values)
    out = cell(1, 0);
else
    out = reshape(sort(values(:)), 1, []);
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


function m = finish_model(ctx, m)
adapter = ModelAdapterManager.getAdapter(ctx.inputs.adapter_matlab);
m = constrainEnzConcs(m);
m = applyKcatConstraints(m);
m = setProtPoolSize(m, [], [], [], adapter);
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


function [ecModel, adapter, fixAdapter] = build_model(ctx)
% The R2/R4-blocked single-route ecTestGEM of flexibilize_edge_cases_ectestgem, concentrations
% filled and not yet constrained.
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
end


function restore_solver(previous)
if ~isempty(previous)
    setRavenSolver(previous);
end
end
