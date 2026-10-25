function results = kcat_merge_edge_cases_ectestgem(ctx)
% MATLAB side of the kcat merge edge-case scenario, mirroring run.py.
%
% Every case builds its own state, runs one call and records what came back and whether the
% call raised. A raise is a value, not a failure. Records are cell arrays of structs;
% lists whose order is the result (the written DLKcat input, the fuzzy matches of a subset)
% are emitted in the order returned.
%
% loadBRENDAdata, fuzzyKcatMatching and the evolutionary distances resolve their files
% through the adapter's own params.path, unlike their Python counterparts, so a *copy* of
% the adapter (ModelAdapter is a value class) is redirected to the scenario's folders.

inputs = ctx.inputs;
adapter = ModelAdapterManager.getAdapter(inputs.adapter_matlab);

results.merge_origin_limit = case_origin_limit();
results.merge_okp_prefix   = case_okp_prefix();
results.assign_zero_row    = case_assign_zero_row(adapter);
results.read_dlkcat        = case_read_dlkcat(adapter, inputs);
results.write_dlkcat       = case_write_dlkcat(adapter, inputs);
results.brenda_outlier     = case_outlier(adapter, inputs);
results.fuzzy              = case_fuzzy(adapter, inputs);
end


% --------------------------------------------------------------------------- %
% helpers
% --------------------------------------------------------------------------- %

function ecModel = build_model(adapter)
model = loadConventionalGEM('modelAdapter', adapter);
ecModel = makeEcModel(model, false, adapter);
ecModel = getECfromGEM(ecModel);
end


function out = num(value)
value = double(value);
if isnan(value)
    out = struct('nan', true, 'value', 0);
else
    out = struct('nan', false, 'value', value);
end
end


function out = sort_records(records, keyFields)
keys = cell(1, numel(records));
for k = 1:numel(records)
    parts = cellfun(@(f) records{k}.(f), keyFields, 'UniformOutput', false);
    keys{k} = strjoin(parts, char(1));
end
[~, order] = sort(keys);
out = records(order);
end


% --------------------------------------------------------------------------- %
% merging kcat lists
% --------------------------------------------------------------------------- %

function list = kcat_list(source, rxns, kcats, wc, origin)
list.source      = source;
list.rxns        = rxns(:);
list.kcats       = kcats(:);
list.wildcardLvl = wc(:);
list.origin      = origin(:);
end


function out = merged_records(list)
n = numel(list.rxns);
records = cell(1, n);
for k = 1:n
    records{k} = struct('reaction', list.rxns{k}, 'kcat', double(list.kcats(k)), ...
        'source', list.kcatSource{k});
end
out = sort_records(records, {'reaction', 'source'});
end


function out = case_origin_limit()
fuzzy = kcat_list('brenda', {'R1'; 'R2'; 'R3'; 'R5'}, [5; 6; 4; 3], zeros(4, 1), [5; 6; 4; 3]);
dlkcat = kcat_list('DLKcat', {'R1'; 'R2'; 'R3'; 'R5'}, [7; 8; 9; 10], nan(4, 1), nan(4, 1));
dlkcat.genes      = repmat({''}, 4, 1);
dlkcat.substrates = repmat({''}, 4, 1);
priority = {'database_top', 'dlkcat', 'database_bottom'};
out.merge_kcats            = merged_records(mergeKcats({fuzzy, dlkcat}, priority));
out.merge_dlkcat_and_fuzzy = merged_records(mergeDLKcatAndFuzzyKcats(dlkcat, fuzzy));
end


function out = case_okp_prefix()
okp = kcat_list('merged', {'R1'; 'R2'}, [3; 4], nan(2, 1), nan(2, 1));
okp.kcatSource = {'OKP-BRENDA'; 'OKP-CataPro'};
dlkcat = kcat_list('DLKcat', {'R1'; 'R2'}, [9; 9], nan(2, 1), nan(2, 1));
priority = {'database_exact', 'catapro', 'dlkcat'};
out.merged = merged_records(mergeKcats({okp, dlkcat}, priority));
end


% --------------------------------------------------------------------------- %
% selecting one kcat per reaction
% --------------------------------------------------------------------------- %

function out = case_assign_zero_row(adapter)
model = build_model(adapter);
list.rxns       = {'R3'; 'R5'};
list.kcats      = [0; 5];
list.kcatSource = {'a'; 'b'};
try
    model = assignKcatValues(model, list, 'criteria', 'max');
    raised = false;
catch
    raised = true;
end
rxns = {'R3', 'R5'};
rows = cell(1, 2);
for k = 1:2
    i = find(strcmp(model.ec.rxns, rxns{k}));
    src = '';
    if ~isempty(model.ec.source{i})
        src = model.ec.source{i};
    end
    rows{k} = struct('reaction', rxns{k}, 'kcat', num(model.ec.kcat(i)), 'source', src);
end
out.raised = raised;
out.rows = rows;
end


% --------------------------------------------------------------------------- %
% DLKcat input and output
% --------------------------------------------------------------------------- %

function out = case_read_dlkcat(adapter, inputs)
names = {'dl_unknown_substrate', 'dl_unknown_rxn_na', 'dl_inf', 'dl_dotted', 'dl_all_none'};
for c = 1:numel(names)
    model = build_model(adapter);
    rows = {};
    try
        list = readDLKcatOutput(model, 'outFile', fullfile(inputs.data_dir, [names{c} '.tsv']), ...
            'modelAdapter', adapter);
        raised = false;
        n = numel(list.rxns);
        records = cell(1, n);
        for k = 1:n
            records{k} = struct('reaction', list.rxns{k}, 'gene', list.genes{k}, ...
                'substrate', list.substrates{k}, 'kcat', num(list.kcats(k)));
        end
        rows = sort_records(records, {'reaction', 'gene'});
    catch
        raised = true;
    end
    out.(names{c}).raised = raised;
    out.(names{c}).rows = rows;
end
end


function out = case_write_dlkcat(adapter, inputs)
labels = {'empty_selection', 'selection_order'};
selections = {{}, {'R5', 'R3'}};
for c = 1:2
    model = build_model(adapter);
    rows = {};
    mask = ismember(model.ec.rxns, selections{c});
    try
        written = writeDLKcatInput(model, 'ecRxns', mask, 'modelAdapter', adapter, ...
            'onlyWithSmiles', false, 'filename', inputs.dlkcat_write_target, 'overwrite', true);
        raised = false;
        n = size(written, 2);
        rows = cell(1, n);
        for k = 1:n
            rows{k} = struct('reaction', written{1, k}, 'gene', written{2, k}, 'substrate', written{3, k});
        end
    catch
        raised = true;
    end
    out.(labels{c}).raised = raised;
    out.(labels{c}).rows = rows;
end
end


% --------------------------------------------------------------------------- %
% BRENDA loading and fuzzy matching
% --------------------------------------------------------------------------- %

function out = case_outlier(adapter, inputs)
outlierAdapter = adapter;
outlierAdapter.params.path = inputs.outlier_root;
[KCATcell, ~] = loadBRENDAdata(outlierAdapter);
records = cell(1, numel(KCATcell{1}));
for k = 1:numel(KCATcell{1})
    records{k} = struct('eccode', KCATcell{1}{k}, 'substrate', KCATcell{2}{k}, ...
        'organism', KCATcell{3}{k}, 'kcat', double(KCATcell{4}(k)));
end
out.kcat_max = sort_records(records, {'eccode', 'substrate', 'organism'});
end


function out = match_records(kcatList)
n = numel(kcatList.rxns);
out = cell(1, n);
for k = 1:n
    ec = '';
    if ~isempty(kcatList.eccodes{k})
        ec = kcatList.eccodes{k};
    end
    origin = -1;
    if ~isnan(kcatList.origin(k))
        origin = double(kcatList.origin(k));
    end
    wc = -1;
    if ~isnan(kcatList.wildcardLvl(k))
        wc = double(kcatList.wildcardLvl(k));
    end
    out{k} = struct('reaction', kcatList.rxns{k}, 'kcat', double(kcatList.kcats(k)), ...
        'eccode', ec, 'origin', origin, 'wildcard_level', wc);
end
end


function out = case_fuzzy(adapter, inputs)
kcatAdapter = adapter;
kcatAdapter.params.path = inputs.scenario_root;

model = build_model(adapter);
model.ec.eccodes{strcmp(model.ec.rxns, 'R3')} = '9.9.9.9';
out.full_wildcard = run_fuzzy(model, find(strcmp(model.ec.rxns, 'R3')), kcatAdapter);

model = build_model(adapter);
model.ec.eccodes{strcmp(model.ec.rxns, 'R2_EXP_1')} = '1.1.1.5;1.1.1.1';
out.tie_between_ec_codes = run_fuzzy(model, find(strcmp(model.ec.rxns, 'R2_EXP_1')), kcatAdapter);

model = build_model(adapter);
ids = {'R3', 'R2_EXP_1', 'R2_EXP_1'};
idx = cellfun(@(r) find(strcmp(model.ec.rxns, r)), ids);
out.subset_order_and_repeat = run_fuzzy(model, idx, kcatAdapter);

model = build_model(adapter);
try
    % [] selects all reactions in MATLAB; an empty index vector selects none, as an empty list does in Python
    list = fuzzyKcatMatching(model, 'ecRxns', zeros(0,1), 'modelAdapter', kcatAdapter);
    raised = false;
    n = numel(list.rxns);
catch
    raised = true;
    n = 0;
end
out.empty_subset.raised = raised;
out.empty_subset.n_matches = n;
end


function out = run_fuzzy(model, ecRxns, kcatAdapter)
matches = {};
try
    list = fuzzyKcatMatching(model, 'ecRxns', ecRxns, 'modelAdapter', kcatAdapter);
    raised = false;
    matches = match_records(list);
catch
    raised = true;
end
out.raised = raised;
out.matches = matches;
end
