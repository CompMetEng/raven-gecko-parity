function results = kcat_io_edge_cases_ectestgem(ctx)
% MATLAB side of the kcat I/O edge-case scenario, mirroring run.py.
%
% Every case builds its own state, runs one call and records what it left behind and whether
% it raised. A raise is a value, not a failure. No case reaches PubChem: every metabolite
% name used is in the vendored smilesDB.tsv of its case.
%
% findMetSmiles, writeOpenKineticsPredictorInput and getStandardKcat resolve their files
% through the adapter's own params.path, so a *copy* of the adapter (ModelAdapter is a value
% class) is redirected to the root of each case.

inputs = ctx.inputs;
adapter = ModelAdapterManager.getAdapter(inputs.adapter_matlab);

results.smiles      = case_smiles(adapter, inputs);
results.okp_input   = case_okp(adapter, inputs);
results.save_source = case_save_source(adapter, inputs);
results.load_source = case_load_source(adapter, inputs);
results.standard    = case_standard(adapter, inputs);
end


% --------------------------------------------------------------------------- %
% helpers
% --------------------------------------------------------------------------- %

function ecModel = build_model(adapter, light, extra)
model = loadConventionalGEM('modelAdapter', adapter);
if extra
    rxnsToAdd.rxns      = {'R6'; 'R7'; 'R8'};
    rxnsToAdd.rxnNames  = {'R6'; 'R7'; 'Biomass Pseudoreaction'};
    rxnsToAdd.equations = {'m1[c] => m2[c]'; 'm1[c] + m2[c] =>'; 'm1[c] => m2[c]'};
    rxnsToAdd.grRules   = {''; ''; ''};
    model = addRxns(model, rxnsToAdd, 3);
end
ecModel = makeEcModel(model, light, adapter);
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


function out = redirected(adapter, root)
out = adapter;
out.params.path = root;
end


% --------------------------------------------------------------------------- %
% findMetSmiles
% --------------------------------------------------------------------------- %

function out = smiles_snapshot(ecModel)
n = numel(ecModel.mets);
out = struct();
[~, order] = sort(ecModel.mets);
for k = order(:)'
    s = '';
    if ~isempty(ecModel.metSmiles{k})
        s = ecModel.metSmiles{k};
    end
    out.(ecModel.mets{k}) = s;
end
end


function out = case_smiles(adapter, inputs)
ecModel = build_model(adapter, false, false);
try
    ecModel = findMetSmiles(ecModel, 'modelAdapter', redirected(adapter, inputs.dup_root), 'verbose', false);
    raised = false;
catch
    raised = true;
end
out.duplicate_cache_name.raised = raised;
out.duplicate_cache_name.smiles = smiles_snapshot(ensure_smiles(ecModel));

ecModel = build_model(adapter, false, false);
ecModel.metNames{strcmp(ecModel.mets, 'm1c')} = '';
try
    ecModel = findMetSmiles(ecModel, 'modelAdapter', redirected(adapter, inputs.empty_root), 'verbose', false);
    raised = false;
catch
    raised = true;
end
out.empty_metabolite_name.raised = raised;
out.empty_metabolite_name.smiles = smiles_snapshot(ensure_smiles(ecModel));
end


function ecModel = ensure_smiles(ecModel)
if ~isfield(ecModel, 'metSmiles')
    ecModel.metSmiles = repmat({''}, numel(ecModel.mets), 1);
end
end


% --------------------------------------------------------------------------- %
% OpenKineticsPredictor input
% --------------------------------------------------------------------------- %

function out = case_okp(adapter, inputs)
okpAdapter = redirected(adapter, inputs.okp_root);
ecModel = build_model(adapter, false, false);
ecModel = findMetSmiles(ecModel, 'modelAdapter', okpAdapter, 'verbose', false);
lines = {};
try
    writeOpenKineticsPredictorInput(ecModel, 'ecRxns', true(numel(ecModel.ec.rxns), 1), ...
        'modelAdapter', okpAdapter, 'onlyWithSmiles', true, 'overwrite', true);
    raised = false;
    csvFile = fullfile(okpAdapter.params.path, 'data', 'OKP.csv');
    fID = fopen(csvFile);
    raw = textscan(fID, '%s %s', 'Delimiter', ',', 'HeaderLines', 1);
    fclose(fID);
    delete(csvFile);
    lines = reshape(sort(strcat(raw{1}, ',', raw{2})), 1, []);
catch
    raised = true;
end
out.all_reactions.raised = raised;
out.all_reactions.lines = lines;
end


% --------------------------------------------------------------------------- %
% kcat source files
% --------------------------------------------------------------------------- %

function out = case_save_source(adapter, inputs)
ecModel = build_model(adapter, false, false);
kcatList.rxns      = {'R3'; 'R5'};
kcatList.kcats     = [25; 30];
kcatList.substrates = {'m1'; 'm2'};
kcatList.eccodes   = {'1.1.2.1'; ''};
kcatList.source    = 'x';
kcatList.detail    = {sprintf('a\tb'); '; x'};
lines = {};
try
    filename = saveKcatSource(ecModel, kcatList, 'edge', 'modelAdapter', adapter, 'folder', inputs.save_dir);
    raised = false;
    text = fileread(filename);
    all = regexp(text, '\r?\n', 'split');
    if ~isempty(all) && isempty(all{end})
        all(end) = [];
    end
    keep = ~startsWith(all, '# generator:') & ~startsWith(all, '# created:');
    lines = all(keep);
catch
    raised = true;
end
out.raised = raised;
out.lines = lines;
end


function out = case_load_source(adapter, inputs)
ecModel = build_model(adapter, false, false);
rows = {};
nInvalid = 0; nUnmatched = 0; nPartial = 0;
try
    loadKcatSource(ecModel, inputs.odd_quote_source, 'modelAdapter', adapter);
    unclosedRaised = false;
catch
    unclosedRaised = true;
end
try
    [kcatList, report] = loadKcatSource(ecModel, inputs.odd_source, 'modelAdapter', adapter);
    raised = false;
    n = numel(kcatList.rxns);
    records = cell(1, n);
    for k = 1:n
        records{k} = struct('reaction', kcatList.rxns{k}, 'kcat', double(kcatList.kcats(k)));
    end
    [~, order] = sort(cellfun(@(r) r.reaction, records, 'UniformOutput', false));
    rows = records(order);
    nInvalid = height(report.invalid);
    nUnmatched = height(report.unmatched);
    nPartial = height(report.partial);
catch
    raised = true;
end
out.raised = raised;
out.rows = rows;
out.n_invalid = nInvalid;
out.n_unmatched = nUnmatched;
out.n_partial = nPartial;
out.unclosed_quote_raised = unclosedRaised;
end


% --------------------------------------------------------------------------- %
% getStandardKcat
% --------------------------------------------------------------------------- %

function out = kcat_rows(ecModel, ids)
records = {};
for k = 1:numel(ids)
    i = find(strcmp(ecModel.ec.rxns, ids{k}));
    if isempty(i)
        continue
    end
    src = '';
    if ~isempty(ecModel.ec.source{i})
        src = ecModel.ec.source{i};
    end
    records{end+1} = struct('reaction', ids{k}, 'kcat', num(ecModel.ec.kcat(i)), 'source', src); %#ok<AGROW>
end
[~, order] = sort(cellfun(@(r) r.reaction, records, 'UniformOutput', false));
out = records(order);
end


function out = case_standard(adapter, inputs)
kcatAdapter = redirected(adapter, inputs.scenario_root);

ecModel = build_model(adapter, false, true);
ecModel.ec.kcat(:) = 10;
before = ecModel.ec.rxns;
try
    ecModel = getStandardKcat(ecModel, 'modelAdapter', kcatAdapter);
    raised = false;
catch
    raised = true;
end
out.extra_reactions.raised = raised;
out.extra_reactions.new_ec_rxns = reshape(sort(setdiff(ecModel.ec.rxns, before)), 1, []);
out.extra_reactions.rows = kcat_rows(ecModel, {'R6', 'R7', 'R8'});

ecModel = build_model(adapter, true, false);
ecModel.ec.kcat(:) = 10;
labels = {'first_call', 'second_call'};
for k = 1:2
    try
        ecModel = getStandardKcat(ecModel, 'modelAdapter', kcatAdapter);
        raised = false;
    catch
        raised = true;
    end
    out.light_called_twice.(labels{k}).raised = raised;
    out.light_called_twice.(labels{k}).enzymes = reshape(sort(ecModel.ec.enzymes), 1, []);
    out.light_called_twice.(labels{k}).ec_rxns = reshape(sort(ecModel.ec.rxns), 1, []);
end
end
