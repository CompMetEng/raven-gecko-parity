function results = merge_compartments_smallyeast(ctx)
% MATLAB side of the compartment-merging scenario.
%
% mergeCompartments groups by metNames. Merged metabolites are compared by id,
% name and ChEBI annotation. smallYeast is merged first, then each toy case, a
% model built from the scenario inputs --- see scenario.yml.
%
% Reactions removed by the merge are reported as the difference between the
% model before and after: RAVEN's deletedRxns covers the deleteRxnsWithOneMet
% mode only, and the reactions that cancelled to nothing are removed separately
% without being reported, so the returned lists are not the same quantity on
% the two sides.

inputs = ctx.inputs;
flag = logical(inputs.delete_single_metabolite_reactions);
results = merge_summary(inputs, readYAMLmodel(inputs.model), flag);
cases = as_list(inputs.toy_cases);
toyResults = cell(1, numel(cases));
for k = 1:numel(cases)
    caseFlag = flag;
    if isfield(cases{k}, 'delete_single_metabolite_reactions')
        caseFlag = logical(cases{k}.delete_single_metabolite_reactions);
    end
    summary = merge_summary(inputs, toy_model(cases{k}), caseFlag);
    toyResults{k} = cell2struct([{cases{k}.name}; struct2cell(summary)], ...
        [{'name'}; fieldnames(summary)], 1);
end
results.toy_cases = toyResults;

end


function model = toy_model(spec)
comps = as_cellstr(spec.comps);
mets = as_list(spec.metabolites);
rxns = as_list(spec.reactions);
nMets = numel(mets);
nRxns = numel(rxns);
model = struct();
model.id = spec.name;
model.comps = comps(:);
model.compNames = comps(:);
model.mets = cellfun(@(m) m.id, mets, 'UniformOutput', false)';
model.metNames = cellfun(@(m) m.name, mets, 'UniformOutput', false)';
model.metComps = cellfun(@(m) find(strcmp(comps, m.comp)), mets)';
model.metFormulas = repmat({''}, nMets, 1);
model.metCharges = zeros(nMets, 1);
model.metMiriams = cell(nMets, 1);
for i = 1:nMets
    if ~isempty(mets{i}.chebi)
        model.metMiriams{i} = struct('name', {{'chebi'}}, 'value', {{mets{i}.chebi}});
    end
end
model.rxns = cellfun(@(r) r.id, rxns, 'UniformOutput', false)';
model.rxnNames = model.rxns;
model.S = sparse(nMets, nRxns);
model.lb = zeros(nRxns, 1);
model.ub = zeros(nRxns, 1);
for j = 1:nRxns
    ids = as_cellstr(rxns{j}.mets);
    coefficients = rxns{j}.coefficients;
    for k = 1:numel(ids)
        model.S(strcmp(model.mets, ids{k}), j) = coefficients(k);
    end
    model.lb(j) = rxns{j}.lb;
    model.ub(j) = rxns{j}.ub;
end
model.rev = double(model.lb < 0);
model.c = zeros(nRxns, 1);
model.b = zeros(nMets, 1);
model.genes = {};
model.grRules = repmat({''}, nRxns, 1);
model.rxnGeneMat = sparse(nRxns, 0);
end


function results = merge_summary(inputs, model, deleteRxnsWithOneMet)
beforeRxns = model.rxns(:);
beforeSpecies = unique(model.metNames(:));
nMetsBefore = numel(model.mets);

[merged, ~, duplicateRxns] = mergeCompartments(model, ...
    'deleteRxnsWithOneMet', deleteRxnsWithOneMet, ...
    'distReverse', true);

afterRxns = merged.rxns(:);

results.n_reactions_before = numel(beforeRxns);
results.n_reactions_after = numel(afterRxns);
results.n_species_before = numel(beforeSpecies);
results.n_metabolites_before = nMetsBefore;
results.n_metabolites_after = numel(merged.mets);
results.compartments = sort_cellstr(merged.comps);
results.reactions = sort_cellstr(afterRxns);
results.removed_reactions = sort_cellstr(setdiff(beforeRxns, afterRxns));
results.deduplicated_reactions = sort_cellstr(flatten_indexed(duplicateRxns));
results.metabolites = sort_cellstr(merged.metNames);
[sortedMets, metOrder] = sort(merged.mets);
metRecords = cell(1, numel(sortedMets));
for k = 1:numel(sortedMets)
    metRecords{k} = struct('id', sortedMets{k}, 'name', merged.metNames{metOrder(k)}, ...
        'chebi', chebi(merged, metOrder(k)));
end
results.metabolite_ids = metRecords;

[sortedRxns, order] = sort(afterRxns);
bounds = cell(1, numel(sortedRxns));
for k = 1:numel(sortedRxns)
    i = order(k);
    bounds{k} = struct( ...
        'reaction', sortedRxns{k}, ...
        'lower_bound', double(merged.lb(i)), ...
        'upper_bound', double(merged.ub(i)));
end
results.bounds = bounds;
results.stoichiometry = stoichiometry(merged);
end


function out = chebi(model, i)
out = '';
if ~isfield(model, 'metMiriams') || isempty(model.metMiriams{i})
    return
end
miriam = model.metMiriams{i};
names = cellstr(miriam.name);
values = cellstr(miriam.value);
out = strjoin(values(strcmp(names, 'chebi')), ';');
end


function records = stoichiometry(model)
% Keyed by metabolite id, matching run.py. char(1)
% joins the sort keys because it sorts below every character an identifier or
% a metabolite name can contain, so the ordering matches Python's tuple
% comparison.
[rowIdx, colIdx, coefficients] = find(model.S);

keys = cell(numel(rowIdx), 1);
for k = 1:numel(rowIdx)
    keys{k} = strjoin({model.rxns{colIdx(k)}, model.mets{rowIdx(k)}}, char(1));
end
[~, order] = sort(keys);

records = cell(1, numel(order));
for k = 1:numel(order)
    i = order(k);
    records{k} = struct( ...
        'reaction', model.rxns{colIdx(i)}, ...
        'metabolite', model.mets{rowIdx(i)}, ...
        'coefficient', double(coefficients(i)));
end
end


function out = flatten_indexed(indexed)
% mergeCompartments documents its third output as the reactions deleted for
% being duplicates, but passes through contractModel's *third* output rather
% than its second: an indexed cell array with one entry per retained reaction,
% holding the ids folded into it separated by semicolons, and empty where
% nothing was. run.py gets a flat list of the ids that went, so the indexed
% form is flattened back to that here.
out = {};
if isempty(indexed)
    return
end
for k = 1:numel(indexed)
    entry = indexed{k};
    if isempty(entry)
        continue
    end
    parts = strsplit(char(entry), ';');
    for j = 1:numel(parts)
        name = strtrim(parts{j});
        if ~isempty(name)
            out{end+1} = name; %#ok<AGROW>
        end
    end
end
end


function out = sort_cellstr(values)
if isempty(values)
    out = {};
    return
end
out = reshape(sort(values(:)), 1, []);
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
