function results = upgrade_gecko3_ectestgem(ctx)
% MATLAB side of the GECKO 3 upgrade scenario, mirroring run.py.
%
% Every case builds its own state in a scratch folder, runs the call and records what it left
% behind and whether it raised. A raise is a value, not a failure. Nothing reaches the network.

inputs = ctx.inputs;
root = char(inputs.scratch);

results.in_memory       = case_in_memory(inputs, root);
results.yaml_round_trip = case_yaml(inputs, root);
results.projects        = case_projects(inputs, root);
results.convert_custom  = case_convert(inputs, root);
end


% --------------------------------------------------------------------------- %
% helpers
% --------------------------------------------------------------------------- %

function project = copy_project(template, target)
if isfolder(target)
    rmdir(target, 's');
end
mkdir(target);
copyfile(char(template), target);
cache = fullfile(target, '__pycache__');
if isfolder(cache)
    rmdir(cache, 's');
end
project = target;
end


function [model, adapter] = build_model(project)
adapter = TomlAdapter(project);
conv = loadConventionalGEM('modelAdapter', adapter);
model = makeEcModel(conv, false, adapter);
end


function idx = protein_rxns(model)
idx = find(startsWith(model.rxns, 'usage_prot_') | strcmp(model.rxns, 'prot_pool_exchange'));
end


function rows = summary(model)
% Pool coefficient, bounds and objective coefficient of the protein usage reactions.
pool = find(strcmp(model.mets, 'prot_pool'));
idx = protein_rxns(model);
[~, order] = sort(model.rxns(idx));
idx = idx(order);
rows = cell(numel(idx), 1);
for k = 1:numel(idx)
    i = idx(k);
    rows{k} = struct('rxn', model.rxns{i}, 'pool_coef', full(model.S(pool, i)), ...
        'lb', model.lb(i), 'ub', model.ub(i), 'c', model.c(i));
end
end


function model = make_reverse(model)
% The GECKO 3 shape: usage and pool exchange reactions reverse, with objective weights.
if ~isfield(model, 'c') || isempty(model.c)
    model.c = zeros(numel(model.rxns), 1);
end
idx = protein_rxns(model);
model.S(:, idx) = -model.S(:, idx);
lb = model.lb(idx);
model.lb(idx) = -model.ub(idx);
model.ub(idx) = -lb;
model.c(strcmp(model.rxns, 'prot_pool_exchange')) = 2;
model.c(strcmp(model.rxns, 'usage_prot_P1')) = -1;
end


function model = make_forward_negative_lb(model)
% A GECKO 4 model whose usage reaction of P1 has a negative lower bound.
if ~isfield(model, 'c') || isempty(model.c)
    model.c = zeros(numel(model.rxns), 1);
end
i = strcmp(model.rxns, 'usage_prot_P1');
model.lb(i) = -5;
model.ub(i) = 1000;
model.c(strcmp(model.rxns, 'prot_pool_exchange')) = 2;
end


function out = sorted_cell(x)
x = sort(x(:));
out = x;
if isempty(out)
    out = {};
end
end


% --------------------------------------------------------------------------- %
% in-memory and file round trip
% --------------------------------------------------------------------------- %

function out = case_in_memory(inputs, root)
project = copy_project(inputs.project_template, fullfile(root, 'memory'));
names = {'reverse', 'forward_negative_lb'};
makers = {@make_reverse, @make_forward_negative_lb};
for k = 1:2
    model = build_model(project);
    model = makers{k}(model);
    before = summary(model);
    [model, changed] = upgradeEcModel(model);
    after = summary(model);
    [model, again] = upgradeEcModel(model);
    out.(names{k}) = struct('before', {before}, 'changed', {sorted_cell(changed)}, ...
        'after', {after}, 'changed_again', {sorted_cell(again)}, 'after_again', {summary(model)});
end
end


function out = case_yaml(inputs, root)
project = copy_project(inputs.project_template, fullfile(root, 'yaml'));
names = {'reverse', 'forward_negative_lb'};
makers = {@make_reverse, @make_forward_negative_lb};
for k = 1:2
    model = build_model(project);
    model = makers{k}(model);
    file = fullfile(project, 'models', [names{k} '.yml']);
    writeYAMLmodel(model, file);
    loaded = readYAMLmodel(file);
    % What the reader returns and which reactions the upgrade turns are not recorded: the
    % two readers differ in whether they turn a reaction with a negative lower bound.
    loaded = upgradeEcModel(loaded);
    out.(names{k}) = struct('after', {summary(loaded)});
end
end


% --------------------------------------------------------------------------- %
% projects
% --------------------------------------------------------------------------- %

function names = listing(project)
today = char(datetime('now', 'Format', 'yyyy-MM-dd'));
d = dir(fullfile(project, '**', '*'));
d = d(~[d.isdir]);
names = cell(numel(d), 1);
for k = 1:numel(d)
    rel = fullfile(d(k).folder, d(k).name);
    rel = rel(numel(project)+2:end);
    names{k} = strrep(strrep(rel, filesep, '/'), today, '<date>');
end
names = sort(names);
if isempty(names)
    names = {};
end
end


function rows = report_rows(report)
rows = cell(height(report), 1);
for k = 1:height(report)
    % MATLAB names the argument 'ecModel', geckopy ec_model.
    rows{k} = struct('step', char(report.step{k}), 'status', char(report.status{k}), ...
        'note', strrep(char(report.note{k}), '(''ecModel'')', '(<ecModel argument>)'));
end
end


function rows = records(file)
% Rows of a kcat source file: rxn, enzymes, kcat, source and note.
rows = {};
if ~isfile(file)
    return
end
lines = strsplit(fileread(file), {'\r\n', '\n'});
header = {};
keys = {};
vals = {};
for k = 1:numel(lines)
    line = lines{k};
    if isempty(line) || line(1) == '#'
        continue
    end
    cells = strsplit(line, '\t', 'CollapseDelimiters', false);
    if isempty(header)
        header = cells;
        continue
    end
    cells(end+1:numel(header)) = {''};
    rec = struct();
    for j = 1:numel(header)
        rec.(header{j}) = cells{j};
    end
    keys{end+1, 1} = strjoin({rec.rxn, rec.enzymes, rec.kcat}, char(1)); %#ok<AGROW>
    vals{end+1, 1} = struct('rxn', rec.rxn, 'enzymes', rec.enzymes, 'kcat', rec.kcat, ...
        'source', rec.source, 'note', rec.note); %#ok<AGROW>
end
[~, order] = sort(keys);
rows = vals(order);
end


function rec = project_case(inputs, root, name, opts)
project = copy_project(inputs.project_template, fullfile(root, name));
if strcmp(opts.dlkcat, 'valid')
    copyfile(char(inputs.dlkcat_valid), fullfile(project, 'data', 'DLKcat.tsv'));
end
[model, ~] = build_model(project);
for k = 1:numel(opts.models)
    m = build_model(project);
    m = opts.makers{k}(m);
    writeYAMLmodel(m, fullfile(project, 'models', [opts.models{k} '.yml']));
end
if opts.no_toml
    delete(fullfile(project, 'model_adapter.toml'));
    copyfile(char(inputs.class_adapter), project);
end
args = {'dryRun', opts.dry_run};
if opts.with_model
    args = [args {'ecModel', model}];
end
raised = false;
rows = {};
report = [];   % evalc assigns it; without this MATLAB reads `report` as a function
try
    evalc('report = upgradeGECKOproject(project, args{:});');
    rows = report_rows(report);
catch
    raised = true;
end
rec.raised = raised;
rec.rows = rows;
if opts.twice && ~raised
    raised2 = false;
    rows2 = {};
    report2 = [];
    try
        evalc('report2 = upgradeGECKOproject(project, ''dryRun'', opts.dry_run, ''ecModel'', model);');
        rows2 = report_rows(report2);
    catch
        raised2 = true; %#ok<NASGU>
    end
    rec.rows_second_run = rows2;
end
rec.files = listing(project);
rec.custom = records(fullfile(project, 'data', 'kcats', 'custom.tsv'));
rec.dlkcat = records(fullfile(project, 'data', 'kcats', 'dlkcat.tsv'));
toml = fullfile(project, 'model_adapter.toml');
first = '';
if isfile(toml)
    lines = strsplit(fileread(toml), {'\r\n', '\n'});
    first = lines{1};
end
rec.toml_first_line = first;
end


function o = options(varargin)
o = struct('models', {{}}, 'makers', {{}}, 'no_toml', false, 'dry_run', true, 'with_model', true, ...
    'twice', false, 'dlkcat', 'valid');
for k = 1:2:numel(varargin)
    o.(varargin{k}) = varargin{k+1};
end
end


function out = case_projects(inputs, root)
out.dry_run_without_model = project_case(inputs, root, 'dry_nomodel', ...
    options('with_model', false, 'models', {'gecko3'}, 'makers', {@make_reverse}));
out.dry_run_with_model = project_case(inputs, root, 'dry_model', ...
    options('models', {'gecko3'}, 'makers', {@make_reverse}));
out.apply_twice = project_case(inputs, root, 'apply', ...
    options('models', {'gecko3'}, 'makers', {@make_reverse}, 'dry_run', false, 'twice', true));
out.forward_model_negative_lb = project_case(inputs, root, 'forward_lb', ...
    options('models', {'forward'}, 'makers', {@make_forward_negative_lb}, 'dry_run', false));
out.stored_dlkcat_foreign_substrates = project_case(inputs, root, 'foreign', ...
    options('dry_run', false, 'dlkcat', 'foreign'));
out.no_toml = project_case(inputs, root, 'no_toml', options('no_toml', true));
end


% --------------------------------------------------------------------------- %
% convertCustomKcats
% --------------------------------------------------------------------------- %

function out = case_convert(inputs, root)
names = {'stock', 'empty_note_and_duplicate'};
keys = {'custom_stock', 'custom_edge'};
for k = 1:2
    project = copy_project(inputs.project_template, fullfile(root, ['convert_' names{k}]));
    [model, adapter] = build_model(project);
    target = fullfile(project, 'data', 'kcats', 'custom.tsv');
    raised = false;
    rows = 0;
    r = [];
    try
        evalc('[~, r] = convertCustomKcats(model, ''customKcats'', char(inputs.(keys{k})), ''outFile'', target, ''modelAdapter'', adapter);');
        rows = r.rows;
    catch
        raised = true;
    end
    out.(names{k}) = struct('raised', raised, 'rows', rows, 'records', {records(target)});
end
end
