function results = kcat_review_stats_ectestgem(ctx)
% MATLAB side of the kcat review statistics scenario, mirroring run.py.
%
% Every case builds its own state, runs the call and records what it returned and whether it
% raised. A raise is a value, not a failure. Nothing reaches the network.

inputs = ctx.inputs;
root = char(inputs.scratch);
loader = loader_stats(inputs, root, inputs.brenda_dir);

results.ec_stats                   = stats_records(loader);
results.review_hand_stats          = review(hand_stats());
results.review_hand_stats_optional = review(hand_stats(), 'include', {'standard-fallback', 'tied-isozymes'});
results.review_hand_stats_top      = review(hand_stats(), 'top', 3);
results.review_loader_stats        = review(loader);
results.lowercase_ec_in_mw         = case_lowercase_mw(inputs, root);
results.measurements               = case_measurements(inputs);
end


% --------------------------------------------------------------------------- %
% helpers
% --------------------------------------------------------------------------- %

function rows = assignments()
% reaction, kcat, source, EC code, leverage (see run.py for what each row is meant to meet)
rows = { ...
    'R_A', 10, 'brenda', '1.1.1.1', 9; ...
    'R_B', 2, 'custom', '1.1.1.1', 8; ...
    'R_C', 6, 'Custom', '1.1.1.1', 7; ...
    'R_D', 5, 'brenda', '4.1.1.1', 6; ...
    'R_F', 40, 'custom', '3.5.1.n3', 5.5; ...
    'R_G1', 77, 'brenda', '', 4; ...
    'R_G2', 77, 'brenda', '', 3.9; ...
    'R_G3', 77, 'brenda', '', 3.8; ...
    'R_G4', 77, 'brenda', '', 3.7; ...
    'R_G5', 77, 'brenda', '', 3.6; ...
    'R_L', 0.001, 'brenda', '', 3; ...
    'R_M', 50000, 'STANDARD', '', 2; ...
    'R_N', 20000, 'standard', '', 1.5; ...
    'R_P_EXP_1', 0.005, 'brenda', '', 1.2; ...
    'R_P_EXP_2', 0.005, 'brenda', '', 1.1};
end


function out = num(value)
value = double(value);
if isnan(value)
    out = struct('nan', true, 'value', 0);
else
    out = struct('nan', false, 'value', value);
end
end


function rows = stats_records(stats)
keys = sort(stats.keys());
rows = cell(numel(keys), 1);
for k = 1:numel(keys)
    s = stats(keys{k});
    values = sort(double(s.values(:)'));
    rows{k} = struct('ec', keys{k}, 'n_kcat', s.nKcat, 'kcat_max', num(s.kcatMax), ...
        'kcat_median', num(s.kcatMedian), 'values', {num2cell(values)}, 'n_sa', s.nSa);
end
end


function stats = hand_stats()
stats = containers.Map('KeyType', 'char', 'ValueType', 'any');
stats('1.1.1.1') = struct('nKcat', 3, 'kcatMax', 10, 'kcatMedian', 2, 'values', [2 6 10], 'nSa', 0);
stats('4.1.1.1') = struct('nKcat', 0, 'kcatMax', NaN, 'kcatMedian', NaN, 'values', [], 'nSa', 1);
stats('3.5.1.n3') = struct('nKcat', 1, 'kcatMax', 40, 'kcatMedian', 4, 'values', [4 40], 'nSa', 0);
end


function out = review(stats, varargin)
a = assignments();
rxns = a(:, 1);
kcats = cell2mat(a(:, 2));
sources = a(:, 3);
ecs = a(:, 4);
lev = cell2mat(a(:, 5));
raised = false;
found = {};
try
    f = reviewKcatAssignment(rxns, kcats, sources, ecs, 'ecStats', stats, 'leverage', lev, varargin{:});
    found = cell(height(f), 1);
    for k = 1:height(f)
        found{k} = struct('rxn', char(f.rxnId{k}), 'ec', char(f.ecCode{k}), 'source', char(f.source{k}), ...
            'kcat', f.kcat(k), 'leverage', f.leverage(k), 'n_isozymes', f.nIsozymes(k), ...
            'other_reactions', f.otherReactions(k), ...
            'checks', {strsplit(char(f.checks{k}), ',')});
    end
catch
    raised = true;
end
out = struct('raised', raised, 'findings', {found});
end


function stats = loader_stats(inputs, root, brendaDir)
project = fullfile(root, 'project');
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
files = dir(fullfile(char(brendaDir), '*.tsv'));
for k = 1:numel(files)
    copyfile(fullfile(files(k).folder, files(k).name), brenda);
end
adapter = TomlAdapter(project);
[kcat, sa] = loadBRENDAdata('modelAdapter', adapter);
stats = ecStatsFromBrenda(kcat, sa);
end


function out = case_lowercase_mw(inputs, root)
% A BRENDA folder whose mw.tsv and sa.tsv hold an EC with a lowercase letter.
out = struct('raised', false, 'stats', {{}});
try
    out.stats = stats_records(loader_stats(inputs, root, inputs.brenda_lower_mw_dir));
catch
    out.raised = true;
end
end


function out = case_measurements(inputs)
names = {'meas_tsv', 'meas_xz', 'meas_both'};
for k = 1:numel(names)
    raised = false;
    rows = {};
    n = 0;
    release = '';
    try
        m = loadBrendaMeasurements('folder', char(inputs.(names{k})));
        n = numel(m.kcat);
        release = m.release;
        keys = cell(n, 1);
        rows = cell(n, 1);
        for i = 1:n
            keys{i} = sprintf('%s%c%020.6f', m.ec{i}, char(1), m.kcat(i));
            rows{i} = struct('ec', m.ec{i}, 'kcat', m.kcat(i), 'log10', m.log10kcat(i));
        end
        [~, order] = sort(keys);
        rows = rows(order);
    catch
        raised = true;
    end
    out.(names{k}) = struct('raised', raised, 'n', n, 'release', release, 'rows', {rows});
end
end
