function results = kcat_assignment_pipeline_ectestgem(ctx)
% MATLAB side of the kcat assignment pipeline scenario, mirroring run.py.
%
% Every case builds its own project folder in a scratch folder, runs the call and records what
% it returned and whether it raised. A raise is a value, not a failure. Nothing reaches the
% network.

inputs = ctx.inputs;
root = char(inputs.scratch);

configs = config_texts();
names = fieldnames(configs);
for k = 1:numel(names)
    results.config.(names{k}) = config_case(inputs, root, names{k}, configs.(names{k}));
end

selectFiles  = struct('custom', 'custom', 'brenda', 'brenda', 'dlkcat', 'dlkcat');
combineFiles = struct('custom', 'custom', 'brenda', 'brenda', 'ec_tree', 'ec_tree', 'homologs', 'homologs');
[select, selectFill, combine] = stage_texts();
results.collect.select_complete = collect_case(inputs, root, 'select_complete', selectFiles, select, false);
results.collect.select_dlkcat_without_rows = collect_case(inputs, root, 'select_empty', ...
    struct('custom', 'custom', 'brenda', 'brenda', 'dlkcat', 'empty_dlkcat'), select, false);
results.collect.select_no_custom_file = collect_case(inputs, root, 'select_nocustom', ...
    struct('brenda', 'brenda', 'dlkcat', 'dlkcat'), select, false);
results.collect.select_isozymes_and_standard = collect_case(inputs, root, 'select_fill', selectFiles, selectFill, false);
results.collect.select_assigned = collect_case(inputs, root, 'select_assign', selectFiles, select, true);
results.collect.combine_complete = collect_case(inputs, root, 'combine_complete', combineFiles, combine, false);
results.collect.combine_homologs_without_rows = collect_case(inputs, root, 'combine_empty', ...
    struct('custom', 'custom', 'brenda', 'brenda', 'ec_tree', 'ec_tree', 'homologs', 'empty_homologs'), combine, false);

results.ec_tree       = tree_case(inputs);
results.ec_tree_model = ec_tree_model_case(inputs, root);
end


% --------------------------------------------------------------------------- %
% texts
% --------------------------------------------------------------------------- %

function [select, selectFill, combine] = stage_texts()
select = sprintf(['[kcat_assignment]\nisozymes = false\nstandard = false\n\n[[kcat_assignment.stage]]\n' ...
    'action = "select"\nsources = ["custom", "database_top", "dlkcat", "database_bottom"]\n' ...
    'top_origin_limit = 4\n']);
selectFill = strrep(strrep(select, 'isozymes = false', 'isozymes = true'), 'standard = false', 'standard = true');
combine = sprintf('[kcat_assignment]\npreset = "ec_tree_homologs"\nisozymes = false\nstandard = false\n');
end


function c = config_texts()
c.default = '';
c.preset_gecko_classic = sprintf('[kcat_assignment]\npreset = "gecko_classic"\n');
c.two_stages = sprintf(['[[kcat_assignment.stage]]\naction = "select"\nsources = ["custom", "database_top"]\n' ...
    'top_origin_limit = 1\n\n[[kcat_assignment.stage]]\naction = "combine"\n' ...
    'sources = ["ec_tree"]\ncriteria = "median"\n']);
c.unknown_key = sprintf('[kcat_assignment]\nbogus = 1\n');
c.stage_sources_string = sprintf('[[kcat_assignment.stage]]\nsources = "custom"\n');
c.preset_empty = sprintf('[kcat_assignment]\npreset = ""\n');
c.criteria_capitalised = sprintf('[[kcat_assignment.stage]]\nsources = ["custom"]\ncriteria = "Median"\n');
c.okp_methods_unknown = sprintf('[kcat_sources.okp]\nmethods = ["Nope"]\n');
c.okp_methods_two = sprintf('[kcat_sources.okp]\nmethods = ["CataPro", "UniKP"]\n');
c.preset_and_stage = sprintf(['[kcat_assignment]\npreset = "okp"\n\n[[kcat_assignment.stage]]\n' ...
    'sources = ["custom"]\n']);
c.brenda_aggregation_capitalised = sprintf('[kcat_sources.brenda]\naggregation = "Median"\n');
c.preset_unknown = sprintf('[kcat_assignment]\npreset = "nonsense"\n');
end


% --------------------------------------------------------------------------- %
% helpers
% --------------------------------------------------------------------------- %

function out = num(value)
value = double(value);
if isnan(value)
    out = struct('nan', true, 'value', 0);
else
    out = struct('nan', false, 'value', value);
end
end


function project = make_project(inputs, root, name, files, toml)
project = fullfile(root, name);
if isfolder(project)
    rmdir(project, 's');
end
mkdir(project);
copyfile(char(inputs.project_template), project);
cache = fullfile(project, '__pycache__');
if isfolder(cache)
    rmdir(cache, 's');
end
brenda = fullfile(project, 'data', 'brenda');
if ~isfolder(brenda)
    mkdir(brenda);
end
f = dir(fullfile(char(inputs.brenda_dir), '*.tsv'));
for k = 1:numel(f)
    copyfile(fullfile(f(k).folder, f(k).name), brenda);
end
if ~isempty(files)
    kcats = fullfile(project, 'data', 'kcats');
    if ~isfolder(kcats)
        mkdir(kcats);
    end
    targets = fieldnames(files);
    for k = 1:numel(targets)
        copyfile(fullfile(char(inputs.kcats_dir), [files.(targets{k}) '.tsv']), ...
            fullfile(kcats, [targets{k} '.tsv']));
    end
end
if ~isempty(toml)
    file = fullfile(project, 'model_adapter.toml');
    text = fileread(file);
    fid = fopen(file, 'w', 'n', 'UTF-8');
    fprintf(fid, '%s\n%s', text, toml);
    fclose(fid);
end
end


function [model, adapter] = build_model(project)
adapter = TomlAdapter(project);
conv = loadConventionalGEM('modelAdapter', adapter);
model = makeEcModel(conv, false, adapter);
model = getECfromGEM(model);
end


function v = plain(value)
% A config value as both sides can state it: double, char or cell of char.
if islogical(value) || isnumeric(value)
    if isscalar(value)
        v = double(value);
    else
        v = num2cell(double(value(:)'));
    end
elseif iscell(value) || (isstring(value) && ~isscalar(value))
    v = cellstr(value(:)')';
    v = v(:)';
else
    v = char(value);
end
end


% --------------------------------------------------------------------------- %
% kcatAssignmentConfig
% --------------------------------------------------------------------------- %

function out = config_case(inputs, root, name, toml)
project = make_project(inputs, root, ['config_' name], [], toml);
out = struct('raised', true, 'preset', '', 'isozymes', false, 'standard', false, 'stages', {{}}, ...
    'sources', {{}}, 'okp_methods', {{}});
try
    adapter = TomlAdapter(project);
    config = kcatAssignmentConfig(adapter.params);
catch
    return
end
stages = cell(numel(config.stages), 1);
for s = 1:numel(config.stages)
    st = config.stages(s);
    stages{s} = struct('action', char(st.action), 'sources', {cellstr(st.sources(:)')}, ...
        'criteria', char(st.criteria), 'top_origin_limit', double(st.topOriginLimit), ...
        'bottom_origin_limit', double(st.bottomOriginLimit), 'wildcard_limit', double(st.wildcardLimit));
end
sources = {};
srcNames = sort(fieldnames(config.sources));
for k = 1:numel(srcNames)
    opts = config.sources.(srcNames{k});
    optNames = sort(fieldnames(opts));
    for j = 1:numel(optNames)
        sources{end+1, 1} = struct('source', srcNames{k}, 'option', optNames{j}, ... %#ok<AGROW>
            'value', {plain(opts.(optNames{j}))});
    end
end
methods = cellstr(config.okpMethods(:)')';
out = struct('raised', false, 'preset', char(config.preset), 'isozymes', logical(config.isozymes), ...
    'standard', logical(config.standard), 'stages', {stages}, 'sources', {sources}, ...
    'okp_methods', {methods(:)'});
end


% --------------------------------------------------------------------------- %
% collectKcats / assignKcats
% --------------------------------------------------------------------------- %

function rows = decision_rows(table)
table = sortrows(table, 'ecRxn');
rows = cell(height(table), 1);
for k = 1:height(table)
    rows{k} = struct('ec_rxn', char(table.ecRxn{k}), 'stage', char(table.stage{k}), ...
        'source', char(table.source{k}), 'kcat', num(table.kcat(k)), ...
        'candidates', table.candidates(k), 'spread', num(table.spreadLog10(k)), ...
        'sd', num(table.sdLog10(k)));
end
end


function out = collect_case(inputs, root, name, files, toml, assign)
project = make_project(inputs, root, name, files, toml);
[model, adapter] = build_model(project);
out = struct('raised', true, 'decisions', {{}}, 'missing', {{}}, 'assigned', {{}});
try
    evalc('[decisions, report] = collectKcats(model, ''modelAdapter'', adapter, ''writeReport'', false);');
catch
    return
end
missing = sort(report.missing(:));
if isempty(missing)
    missing = {};
end
out = struct('raised', false, 'decisions', {decision_rows(decisions)}, 'missing', {missing}, ...
    'assigned', {{}});
if assign
    assignRaised = false;
    assigned = {};
    try
        evalc('[model, ~, ~] = assignKcats(model, ''modelAdapter'', adapter, ''writeReport'', false);');
        [rxns, order] = sort(model.ec.rxns);
        kcats = model.ec.kcat(order);
        sources = model.ec.source(order);
        assigned = cell(numel(rxns), 1);
        for k = 1:numel(rxns)
            assigned{k} = struct('ec_rxn', rxns{k}, 'kcat', num(kcats(k)), 'source', char(sources{k}));
        end
    catch
        assignRaised = true;
    end
    out.assign_raised = assignRaised;
    out.assigned = assigned;
end
end


% --------------------------------------------------------------------------- %
% EC tree
% --------------------------------------------------------------------------- %

function q = queries()
q = {'1.1.1.1', '1.1.1.1;2.7.1.2', '', '9.9.9.9', '1.1.2.1;1.1.2.1', '3.5.1.n3', 'x.y', ...
    '2.7.1.1;2.7.1.2', '1.1.1.1;1.1.1.2', '1.1.1.2;1.1.1.9', '2.7.1.2;4.1.1.1', ' 1.1.1.1 ; 2.7.1.1 '};
end


function rows = estimates(measurements, k)
q = queries();
[v, n, depth, ec] = ecTreeEstimate(q(:), measurements, k);
rows = cell(numel(q), 1);
for i = 1:numel(q)
    rows{i} = struct('query', q{i}, 'log10', num(v(i)), 'n', n(i), 'depth', depth(i), 'ec', ec{i});
end
end


function out = tree_case(inputs)
measurements = loadBrendaMeasurements('folder', char(inputs.brenda_dir));
labels = {'k5', 'k1', 'k0', 'k_negative', 'k_nan'};
ks = [5 1 0 -1 NaN];
for i = 1:numel(labels)
    raised = false;
    rows = {};
    try
        rows = estimates(measurements, ks(i));
    catch
        raised = true;
    end
    out.(labels{i}) = struct('raised', raised, 'estimates', {rows});
end
nanTable = struct('ec', {{'1.1.1.1'; '1.1.1.1'; '2.7.1.1'; '2.7.1.1'; '2.7.1.1'}}, ...
    'log10kcat', [1; NaN; 2; 3; 4]);
raised = false;
rows = {};
try
    [v, n, depth, ec] = ecTreeEstimate({'1.1.1.1'; '2.7.1.1'; '9.9.9.9'}, nanTable, 5);
    q = {'1.1.1.1', '2.7.1.1', '9.9.9.9'};
    rows = cell(3, 1);
    for i = 1:3
        rows{i} = struct('query', q{i}, 'log10', num(v(i)), 'n', n(i), 'depth', depth(i), 'ec', ec{i});
    end
catch
    raised = true;
end
out.nan_in_table = struct('raised', raised, 'estimates', {rows});
end


function rows = tree_rows(kcatList)
n = numel(kcatList.rxns);
keys = cell(n, 1);
rows = cell(n, 1);
for i = 1:n
    keys{i} = sprintf('%s%c%030.12f', kcatList.rxns{i}, char(1), kcatList.kcats(i));
    rows{i} = struct('rxn', kcatList.rxns{i}, 'eccode', kcatList.eccodes{i}, ...
        'kcat', num(kcatList.kcats(i)), 'detail', kcatList.detail{i});
end
[~, order] = sort(keys);
rows = rows(order);
end


function out = ec_tree_model_case(inputs, root)
project = make_project(inputs, root, 'ec_tree_model', [], '');
[model, adapter] = build_model(project);
ids = sort(model.ec.rxns);
[~, i0] = ismember(ids{1}, model.ec.rxns);
[~, i1] = ismember(ids{2}, model.ec.rxns);
[~, i2] = ismember(ids{3}, model.ec.rxns);
variants = {'all', {}; 'k1', {'k', 1}; 'subset', {'ecRxns', [i0; i1; i2]}; ...
    'subset_with_repeat', {'ecRxns', [i0; i1; i1]}; 'empty_subset', {'ecRxns', []}};
out.n_ec_rxns = numel(ids);
for v = 1:size(variants, 1)
    raised = false;
    rows = {};
    kcatList = [];
    try
        evalc('kcatList = kcatsFromECTree(model, ''save'', false, ''modelAdapter'', adapter, variants{v, 2}{:});');
        rows = tree_rows(kcatList);
    catch
        raised = true;
    end
    out.(variants{v, 1}) = struct('raised', raised, 'rows', {rows});
end
end
