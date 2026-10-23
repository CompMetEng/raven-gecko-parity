function results = kcat_apply_edge_cases_ectestgem(ctx)
% MATLAB side of the kcat-apply edge-case scenario, mirroring run.py.
%
% Every case builds its own state, runs one call and records the kcat rows, the enzyme
% coefficients of the S matrix and whether the call raised. A raise is a value, not a
% failure. Records are cell arrays of structs (a 1x1 struct array would encode as a bare
% object), sorted on the keys joined with char(1).

inputs = ctx.inputs;
adapter = ModelAdapterManager.getAdapter(inputs.adapter_matlab);

results.light_subset             = case_light_subset();
results.light_empty_isozyme_row  = case_light_empty_row();
results.light_update_ids         = case_light_ids();
results.full_zero_reapply        = case_full_zero_reapply(adapter);
results.full_update_list         = case_full_update_list(adapter);
results.set_kcat                 = case_set_kcat(adapter);
results.isozymes                 = case_isozymes(adapter);
results.custom                   = case_custom(adapter, inputs.custom_dir);
end


% --------------------------------------------------------------------------- %
% helpers
% --------------------------------------------------------------------------- %

function ecModel = build_model(adapter, light)
model = loadConventionalGEM('modelAdapter', adapter);
ecModel = makeEcModel(model, light, adapter);
end


function out = prot_coefficients(model)
protIdx = find(startsWith(model.mets, 'prot_'));
% Coefficients of the metabolic reactions only: the usage reactions and the pool exchange
% carry fixed -1 / +1 entries that do not depend on any kcat.
records = {};
for p = 1:numel(protIdx)
    row = protIdx(p);
    [~, rxnIdx, vals] = find(model.S(row, :));
    for k = 1:numel(rxnIdx)
        if startsWith(model.rxns{rxnIdx(k)}, 'usage_prot_') || strcmp(model.rxns{rxnIdx(k)}, 'prot_pool_exchange')
            continue
        end
        records{end+1} = struct('reaction', model.rxns{rxnIdx(k)}, ...
            'metabolite', model.mets{row}, 'coefficient', full(double(vals(k)))); %#ok<AGROW>
    end
end
keys = cellfun(@(r) [r.reaction char(1) r.metabolite], records, 'UniformOutput', false);
[~, order] = sort(keys);
out = records(order);
end


function out = num(value)
value = double(value);
if isnan(value)
    out = struct('nan', true, 'value', 0);
else
    out = struct('nan', false, 'value', value);
end
end


function out = kcat_rows(ecModel)
n = numel(ecModel.ec.rxns);
records = cell(1, n);
for i = 1:n
    src = '';
    if ~isempty(ecModel.ec.source{i})
        src = ecModel.ec.source{i};
    end
    nt = '';
    if ~isempty(ecModel.ec.notes{i})
        nt = ecModel.ec.notes{i};
    end
    records{i} = struct('reaction', ecModel.ec.rxns{i}, 'kcat', num(ecModel.ec.kcat(i)), ...
        'source', src, 'notes', nt);
end
[~, order] = sort(ecModel.ec.rxns);
out = records(order);
end


function model = set_kcats(model, ids, values)
for k = 1:numel(ids)
    model.ec.kcat(strcmp(model.ec.rxns, ids{k})) = values(k);
end
end


% --------------------------------------------------------------------------- %
% hand-built light models
% --------------------------------------------------------------------------- %

function model = light_model(kcat, enzRows, mw)
clear model
n = numel(kcat);
model.rxns = {'R1'};
model.mets = {'prot_pool'};
model.S = sparse(1, 1);
model.lb = 0;
model.ub = 1000;
model.ec.geckoLight = true;
model.ec.rxns = arrayfun(@(i) sprintf('%03d_R1', i), (1:n)', 'UniformOutput', false);
model.ec.kcat = kcat(:);
model.ec.mw = mw(:);
model.ec.rxnEnzMat = sparse(enzRows);
end


function out = case_light_subset()
% MATLAB names the cobra reaction R1, so all isozymes compete; Python names the ec rows
% it changed.
labels = {'changed_isozyme_only', 'both_isozymes'};
for k = 1:2
    model = light_model([200; 10], [1 0; 0 1], [10000; 5000]);
    model = applyKcatConstraints(model);
    model.ec.kcat(2) = 20;
    model = applyKcatConstraints(model, {'R1'});
    out.(labels{k}) = prot_coefficients(model);
end
end


function out = case_light_empty_row()
% 002_R1 has a kcat but no enzyme in the coupling matrix.
model = light_model([50; 30], [1 0; 0 0], [10000; 5000]);
model = applyKcatConstraints(model);
out.coefficients = prot_coefficients(model);
end


function out = case_light_ids()
ids = {'cobra_id', 'ec_id'};
values = {{'R1'}, {'001_R1'}};
for k = 1:2
    model = light_model([200; 10], [1 0; 0 1], [10000; 5000]);
    try
        applyKcatConstraints(model, values{k});
        raised = false;
    catch
        raised = true;
    end
    out.(ids{k}).raised = raised;
end
end


% --------------------------------------------------------------------------- %
% applyKcatConstraints on the full model
% --------------------------------------------------------------------------- %

function out = case_full_zero_reapply(adapter)
model = build_model(adapter, false);
model.ec.kcat(:) = 10;
model = applyKcatConstraints(model);
out.after_first_apply = prot_coefficients(model);
model.ec.kcat(:) = 0;
model = applyKcatConstraints(model);
out.after_zero_reapply = prot_coefficients(model);
end


function out = case_full_update_list(adapter)
model = build_model(adapter, false);
model.ec.kcat(:) = 10;
model = applyKcatConstraints(model);
model.ec.kcat(:) = 20;
try
    model = applyKcatConstraints(model, 'updateRxns', {});
    raised = false;
catch
    raised = true;
end
out.empty_list.raised = raised;
out.empty_list.coefficients = prot_coefficients(model);

model = build_model(adapter, false);
model.ec.kcat(:) = 10;
try
    model = applyKcatConstraints(model, {'NOSUCH'});
    raised = false;
catch
    raised = true;
end
out.unknown_id.raised = raised;
out.unknown_id.coefficients = prot_coefficients(model);
end


% --------------------------------------------------------------------------- %
% setKcatForReactions
% --------------------------------------------------------------------------- %

function out = case_set_kcat(adapter)
labels = {'default_apply', 'lowercase_id', 'unknown_id'};
ids = {'R3', 'r3', 'NOSUCH'};
for k = 1:3
    model = build_model(adapter, false);
    try
        model = setKcatForReactions(model, ids{k}, 100);
        raised = false;
    catch
        raised = true;
    end
    out.(labels{k}).raised = raised;
    out.(labels{k}).rows = kcat_rows(model);
    out.(labels{k}).coefficients = prot_coefficients(model);
end
model = build_model(adapter, true);
try
    model = setKcatForReactions(model, 'R3', 100);
    raised = false;
catch
    raised = true;
end
out.light_base_id.raised = raised;
out.light_base_id.rows = kcat_rows(model);
end


% --------------------------------------------------------------------------- %
% getKcatAcrossIsozymes
% --------------------------------------------------------------------------- %

function out = case_isozymes(adapter)
model = build_model(adapter, false);
model = set_kcats(model, {'R2_EXP_1', 'R2_EXP_2', 'R2_REV_EXP_1', 'R2_REV_EXP_2', 'R3'}, [NaN 5 0 NaN 8]);
try
    model = getKcatAcrossIsozymes(model);
    raised = false;
catch
    raised = true;
end
out.nan_siblings.raised = raised;
out.nan_siblings.rows = kcat_rows(model);

model = build_model(adapter, false);
model = set_kcats(model, {'R2_EXP_1', 'R2_EXP_2', 'R3'}, [4 0 8]);
try
    model = getKcatAcrossIsozymes(model);
    raised = false;
catch
    raised = true;
end
out.default_apply.raised = raised;
out.default_apply.rows = kcat_rows(model);
out.default_apply.coefficients = prot_coefficients(model);
end


% --------------------------------------------------------------------------- %
% applyCustomKcats
% --------------------------------------------------------------------------- %

function out = case_custom(adapter, folder)
names = {'ck_empty_notes', 'ck_blank_kcat', 'ck_lowercase', 'ck_protein_only', ...
    'ck_quoted_list', 'ck_partial_complex'};
for k = 1:numel(names)
    model = build_model(adapter, false);
    try
        model = applyCustomKcats(model, 'customKcats', fullfile(folder, [names{k} '.tsv']), ...
            'modelAdapter', adapter);
        raised = false;
    catch
        raised = true;
    end
    out.(names{k}).raised = raised;
    out.(names{k}).rows = kcat_rows(model);
end
labels = {'light_model', 'light_protein_only'};
files = {'ck_empty_notes', 'ck_protein_only'};
for k = 1:2
    model = build_model(adapter, true);
    try
        model = applyCustomKcats(model, 'customKcats', fullfile(folder, [files{k} '.tsv']), ...
            'modelAdapter', adapter);
        raised = false;
    catch
        raised = true;
    end
    out.(labels{k}).raised = raised;
    out.(labels{k}).rows = kcat_rows(model);
end
end
