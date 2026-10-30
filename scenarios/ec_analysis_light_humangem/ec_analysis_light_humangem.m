function results = ec_analysis_light_humangem(ctx)
% MATLAB side of the light Human-GEM ecModel analysis scenario; mirrors run.py. Solving,
% mapping and subsetting on the tutorial's gecko-light ecModel. See
% scenario.yml.

previousSolver = getpref('RAVEN', 'solver', '');
restoreSolver = onCleanup(@() restore_solver(previousSolver));
setRavenSolver(char(ctx.inputs.matlab_solver));

inputs = ctx.inputs;
tutorialAdapter = ModelAdapterManager.getAdapter(inputs.adapter_matlab);
conv = loadConventionalGEM('modelAdapter', tutorialAdapter);
bioId = tutorialAdapter.getParameters().bioRxn;

% The tutorial ecModel has 52 ec reactions with no enzyme written as an empty `enzymes: !!omap`
% entry; both readers load the file as it is (see ec_model_load_light_humangem). Both sides here read a
% copy without them; loadEcModel joins its filename onto <adapter path>/models/, so the copy
% goes to a temporary folder that a copy of the adapter points to (ModelAdapter is a value class).
workdir = tempname;
mkdir(fullfile(workdir, 'models'));
cleanup = onCleanup(@() rmdir(workdir, 's'));
text = fileread(fullfile(tutorialAdapter.getParameters().path, 'models', char(inputs.model_file)));
text = regexprep(text, '(?m)^ *- enzymes: !!omap\n(?= {0,4}\S)', '');
fid = fopen(fullfile(workdir, 'models', 'repaired.yml'), 'w', 'n', 'UTF-8');
fprintf(fid, '%s', text);
fclose(fid);
adapter = tutorialAdapter;
adapter.params.path = workdir;
inputs.model_file = 'repaired.yml';
model = fresh(inputs, adapter, bioId);
results = struct();

sol = solveLP(model);
results.growth = run_case(@() struct('objective', num(sol.f), ...
    'pool', num(sol.x(strcmp(model.rxns, 'prot_pool_exchange')))));
results.min_pool = run_case(@() min_pool_case(inputs, adapter, bioId));
results.mapping = run_case(@() mapping_case(model, conv));
results.subset = run_case(@() subset_case(inputs, adapter, bioId));
results.pool_sizes = run_case(@() pool_sizes_case(inputs, adapter, bioId));

end


function model = fresh(inputs, adapter, bioId)
model = loadEcModel(char(inputs.model_file), adapter);
model.c = double(strcmp(model.rxns, bioId));
end


function out = min_pool_case(inputs, adapter, bioId)
m = fresh(inputs, adapter, bioId);
optimum = solveLP(m);
bioIdx = strcmp(m.rxns, bioId);
m.lb(bioIdx) = 0.99 * optimum.f;
poolIdx = strcmp(m.rxns, 'prot_pool_exchange');
m.c = zeros(numel(m.rxns), 1);
m.c(poolIdx) = -1;
s = solveLP(m);
out = struct('pool_usage', num(-s.f), 'growth_lower_bound', num(m.lb(bioIdx)));
end


function out = mapping_case(model, conv)
% A fixed flux vector, (i * 7919 mod 1013) / 1013 - 0.5 for the i-th reaction (from 1), given to
% both sides: the LP solution differs between the two (see scenario.yml).
idx = (1:numel(model.rxns)).';
vector = mod(idx * 7919, 1013) / 1013 - 0.5;
[mappedFlux, enzUsageFlux, usageEnz] = mapRxnsToConv(model, conv, vector);
n = numel(conv.rxns);
mapped = cell(1, n);
for k = 1:n
    mapped{k} = struct('reaction', conv.rxns{k}, 'flux', num(mappedFlux(k)));
end
m = numel(usageEnz);
usage = cell(1, m);
for k = 1:m
    usage{k} = struct('label', usageEnz{k}, 'flux', num(enzUsageFlux(k)));
end
[~, order] = sort(usageEnz);
out = struct('mapped', {mapped}, 'usage', {usage(order)});
end


function out = subset_case(inputs, adapter, bioId)
small = readYAMLmodel(char(inputs.subset_file));
smallGrowth = solveLP(small).f;
ecModel = fresh(inputs, adapter, bioId);
sub = getSubsetEcModel(ecModel, small);
sub.c = double(strcmp(sub.rxns, bioId));
out = struct('small_growth', num(smallGrowth), 'ec_growth', num(solveLP(ecModel).f), ...
    'subset_growth', num(solveLP(sub).f));
out.counts = struct('reactions', numel(sub.rxns), 'metabolites', numel(sub.mets), ...
    'genes', numel(sub.genes), 'ec_rxns', numel(sub.ec.rxns), 'enzymes', numel(sub.ec.enzymes));
out.reactions = reshape(sort(sub.rxns), 1, []);
out.ec_rxns = reshape(sub.ec.rxns, 1, []);
out.enzymes = reshape(sub.ec.enzymes, 1, []);
end


function out = pool_sizes_case(inputs, adapter, bioId)
m = fresh(inputs, adapter, bioId);
m = setProtPoolSize(m, [], [], [], adapter);
defaultUb = m.ub(strcmp(m.rxns, 'prot_pool_exchange'));
params = adapter.getParameters();
m = setProtPoolSize(m, 'Ptot', params.Ptot * double(inputs.pool_scale), 'modelAdapter', adapter);
halvedUb = m.ub(strcmp(m.rxns, 'prot_pool_exchange'));
s = solveLP(m);
out = struct('default_ub', num(defaultUb), 'halved_ub', num(halvedUb), 'halved_growth', num(s.f));
end


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
catch err
    fprintf('case raised: %s\n', err.message);
    out = struct('raised', true, 'value', struct());
    return
end
out = struct('raised', false, 'value', {value});
end


function restore_solver(previous)
if ~isempty(previous)
    setRavenSolver(previous);
end
end
