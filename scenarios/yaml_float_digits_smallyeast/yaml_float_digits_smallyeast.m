function results = yaml_float_digits_smallyeast(ctx)
% MATLAB side of the YAML float-digits scenario.
%
% Sets bounds on three reactions, writes the model and returns the lines of the
% file that hold a lower_bound or upper_bound. The file goes to a temporary
% folder: what the writer prints is the point, not an artefact to keep.

inputs = ctx.inputs;
model = readYAMLmodel(inputs.model);

ids = fieldnames(inputs.bounds);
for k = 1:numel(ids)
    i = find(strcmp(model.rxns, ids{k}), 1);
    model.lb(i) = double(inputs.bounds.(ids{k}).lb);
    model.ub(i) = double(inputs.bounds.(ids{k}).ub);
end

workdir = tempname;
mkdir(workdir);
cleanup = onCleanup(@() rmdir(workdir, 's'));
target = fullfile(workdir, 'model.yml');
writeYAMLmodel(model, 'fileName', target);

text = strrep(fileread(target), sprintf('\r\n'), sprintf('\n'));
lines = strsplit(text, sprintf('\n'), 'CollapseDelimiters', false);
keep = contains(lines, 'lower_bound') | contains(lines, 'upper_bound');
lines = lines(keep);

results.n_bound_lines = numel(lines);
results.bound_lines = reshape(lines, 1, []);

end
