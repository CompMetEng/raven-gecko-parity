function results = export_for_git_toy(ctx)
% MATLAB side of the exportForGit scenario.
%
% A four-reaction model is exported twice: as txt only, to compare the table, and with the
% default formats, to compare which files are written. Equations are recorded apart from the
% rest of each row, because the two sides spell them differently.

reactions = as_list(ctx.inputs.reactions);
n = numel(reactions);
mets = {};
for i = 1:n
    mets = [mets, as_cellstr(reactions{i}.substrates), as_cellstr(reactions{i}.products)]; %#ok<AGROW>
end
mets = unique(mets, 'stable');
S = zeros(numel(mets), n);
model = struct();
model.id = 'exportToy';
model.name = 'exportToy';
model.version = '1.0.0';
model.date = '2026-10-22';
model.rxns = cell(n, 1);
model.rxnNames = cell(n, 1);
model.lb = zeros(n, 1);
model.ub = zeros(n, 1);
model.c = zeros(n, 1);
model.grRules = cell(n, 1);
for i = 1:n
    r = reactions{i};
    model.rxns{i} = char(r.id);
    model.rxnNames{i} = char(r.id);
    model.lb(i) = double(r.lb);
    model.ub(i) = double(r.ub);
    model.c(i) = double(r.objective);
    model.grRules{i} = char(r.rule);
    for m = as_cellstr(r.substrates)
        S(strcmp(mets, m{1}), i) = -1;
    end
    for m = as_cellstr(r.products)
        S(strcmp(mets, m{1}), i) = 1;
    end
end
model.mets = strcat(mets(:), '_c');
model.metNames = mets(:);
model.metFormulas = repmat({'C'}, numel(mets), 1);
model.S = sparse(S);
model.rev = double(model.lb < 0);
model.b = zeros(numel(mets), 1);
model.comps = {'c'};
model.compNames = {'cytosol'};
model.metComps = ones(numel(mets), 1);
genes = {};
for i = 1:n
    tokens = regexp(model.grRules{i}, '[A-Za-z0-9_.]+', 'match');
    genes = [genes, tokens(~ismember(lower(tokens), {'and', 'or'}))]; %#ok<AGROW>
end
model.genes = unique(genes, 'stable')';
model.rxnGeneMat = sparse(n, numel(model.genes));
for i = 1:n
    tokens = regexp(model.grRules{i}, '[A-Za-z0-9_.]+', 'match');
    tokens = tokens(~ismember(lower(tokens), {'and', 'or'}));
    [~, col] = ismember(tokens, model.genes);
    model.rxnGeneMat(i, col) = 1;
end

workdir = tempname;
mkdir(workdir);
cleanup = onCleanup(@() rmdir(workdir, 's'));

txtDir = fullfile(workdir, 'txt');
evalc('exportForGit(model, ''path'', txtDir, ''formats'', {''txt''})');
lines = splitlines(strtrim(fileread(fullfile(txtDir, 'model', 'txt', 'model.txt'))));
dependencies = splitlines(strtrim(fileread(fullfile(txtDir, 'model', 'dependencies.txt'))));

allDir = fullfile(workdir, 'all');
evalc('exportForGit(model, ''path'', allDir)');
found = dir(fullfile(allDir, 'model', '**', '*'));
found = found(~[found.isdir]);
written = cell(1, numel(found));
for k = 1:numel(found)
    [~, parent] = fileparts(found(k).folder);
    written{k} = [parent '/' found(k).name];
end

rows = cell(1, numel(lines) - 1);
for k = 2:numel(lines)
    fields = strsplit(lines{k}, sprintf('\t'), 'CollapseDelimiters', false);
    rows{k - 1} = struct('reaction', fields{1}, 'equation', fields{2}, 'rule', fields{3}, ...
        'lower_bound', fields{4}, 'upper_bound', fields{5}, 'objective', fields{6});
end
[~, order] = sort(cellfun(@(r) r.reaction, rows, 'UniformOutput', false));
rows = rows(order);

names = cellfun(@(l) strtok(l, sprintf('\t')), dependencies, 'UniformOutput', false);
results.header = strsplit(lines{1}, sprintf('\t'), 'CollapseDelimiters', false);
results.rows = rows;
results.libsbml_recorded = any(strcmp(names, 'libSBML'));
results.files_written = reshape(sort(written), 1, []);
results.flag_checks = flag_checks(model, as_list(ctx.inputs.flag_cases), workdir);
end


function checks = flag_checks(model, flagCases, workdir)
% Whether a call with one flag value is refused as RAVEN:badInput; text is refused, a logical is not.
names = struct('main_branch_flag', 'mainBranchFlag', 'sub_dirs', 'subDirs');
checks = struct();
for k = 1:numel(flagCases)
    c = flagCases{k};
    target = fullfile(workdir, char(c.name));
    mkdir(target);
    try
        evalc('exportForGit(model, ''path'', target, ''formats'', {''txt''}, names.(char(c.flag)), c.value)');
        checks.(char(c.name)) = struct('rejected', false);
    catch ME
        if ~strcmp(ME.identifier, 'RAVEN:badInput')
            rethrow(ME);
        end
        checks.(char(c.name)) = struct('rejected', true);
    end
end
end


function out = as_cellstr(value)
if isempty(value)
    out = {};
elseif ischar(value)
    out = {value};
else
    out = reshape(cellfun(@char, value, 'UniformOutput', false), 1, []);
end
end


function out = as_list(value)
% jsondecode gives a struct array when every element has the same fields and a cell otherwise.
if iscell(value)
    out = value(:)';
else
    out = num2cell(value(:)');
end
end
