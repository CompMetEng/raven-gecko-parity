function results = relax_enz_concs_ectestgem(ctx)
% MATLAB side of the relaxEnzConcs scenario; mirrors run.py. relaxEnzConcs with perComplex on
% ecTestGEM with a complex written into ec.rxnEnzMat. A case that raises is a result, not a
% failure.

previousSolver = getpref('RAVEN', 'solver', '');
restoreSolver = onCleanup(@() restore_solver(previousSolver));
setRavenSolver(char(ctx.inputs.matlab_solver));

inputs = ctx.inputs;
results = struct();
labels = {'measured', 'median', 'high'};
statistics = {[], 'median', []};
targets = [double(inputs.exp_growth), double(inputs.exp_growth), double(inputs.exp_growth_high)];
for k = 1:numel(labels)
    for perComplex = [false true]
        if perComplex
            suffix = 'per_complex';
        else
            suffix = 'default';
        end
        name = ['relax_' labels{k} '_' suffix];
        results.(name) = run_case(@() relax_case(ctx, statistics{k}, targets(k), perComplex));
    end
end
end


function out = relax_case(ctx, statistic, expGrowth, perComplex)
% An empty statistic leaves the measured concentrations of the complex as they are.
inputs = ctx.inputs;
[m, ~, fixAdapter] = build_model(ctx);
m = set_layout(m, inputs.layouts.binding);
if ~isempty(statistic)
    m = levelComplexConcs(m, 'statistic', statistic);
end
m = finish_model(ctx, m);
start = solveLP(m);
[mRelax, flexEnz] = relaxEnzConcs(m, 'expGrowth', expGrowth, 'bioRxn', 'R5', ...
    'perComplex', perComplex, 'modelAdapter', fixAdapter);
ids = reshape(flexEnz.uniprotIDs, 1, []);
n = numel(ids);
records = cell(1, n);
for i = 1:n
    records{i} = struct('protein', ids{i}, 'old_conc', num(flexEnz.oldConcs(i)), ...
        'flex_conc', num(flexEnz.flexConcs(i)), 'ratio_incr', num(flexEnz.ratioIncr(i)));
end
if n > 0
    [~, order] = sort(cellfun(@(r) r.protein, records, 'UniformOutput', false));
    records = records(order);
end
final = solveLP(mRelax);
optimal = final.stat == 1;
growth = 0;
if optimal
    growth = final.f;
end
out = struct('concs', {concs_list(m)}, 'start_growth', num(start.f), 'relaxed', {records}, ...
    'relaxed_set', {sorted_row(ids)}, 'final_is_optimal', logical(optimal), ...
    'final_growth', num(growth), 'usage_ub', {usage_ub_list(mRelax)}, ...
    'pool_ub', num(mRelax.ub(strcmp(mRelax.rxns, 'prot_pool_exchange'))));
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
