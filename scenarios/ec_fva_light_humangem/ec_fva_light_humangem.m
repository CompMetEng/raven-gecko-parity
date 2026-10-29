function results = ec_fva_light_humangem(ctx)
% MATLAB side of the light Human-GEM flux variability scenario; mirrors run.py. ecFVA on the
% tutorial's gecko-light ecModel with growth fixed at its optimum. See scenario.yml.
%
% The tutorial ecModel has 52 ec reactions with no enzyme written as an empty `enzymes: !!omap`
% entry, which geckopy's reader rejects (see ec_model_load_light_humangem). Both sides read a
% copy without them; loadEcModel joins its filename onto <adapter path>/models/, so the copy
% goes to a temporary folder that a copy of the adapter points to (ModelAdapter is a value
% class).

previousSolver = getpref('RAVEN', 'solver', '');
restoreSolver = onCleanup(@() restore_solver(previousSolver));
setRavenSolver(char(ctx.inputs.matlab_solver));

inputs = ctx.inputs;
tutorialAdapter = ModelAdapterManager.getAdapter(inputs.adapter_matlab);
conv = loadConventionalGEM('modelAdapter', tutorialAdapter);
bioId = tutorialAdapter.getParameters().bioRxn;

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

m = loadEcModel('repaired.yml', adapter);
m.c = double(strcmp(m.rxns, bioId));
optimum = solveLP(m);
bioIdx = strcmp(m.rxns, bioId);
m.lb(bioIdx) = optimum.f;
m.ub(bioIdx) = optimum.f;
[minFlux, maxFlux] = ecFVA(m, conv);

stride = double(inputs.fva_stride);
idx = 1:stride:numel(conv.rxns);
rows = cell(1, numel(idx));
for k = 1:numel(idx)
    rows{k} = struct('reaction', conv.rxns{idx(k)}, 'min_flux', num(minFlux(idx(k))), ...
        'max_flux', num(maxFlux(idx(k))));
end
results = struct('optimum', num(optimum.f), 'n_reactions', double(numel(conv.rxns)), 'sample', {rows});

end


function out = num(x)
% A number as a record, so that NaN is a value both sides can state.
x = double(x);
if isnan(x)
    out = struct('is_nan', true, 'value', 0);
else
    out = struct('is_nan', false, 'value', x);
end
end


function restore_solver(previous)
if ~isempty(previous)
    setRavenSolver(previous);
end
end
