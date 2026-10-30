function results = ec_model_build_edge_cases_ectestgem(ctx)
% MATLAB side of the makeEcModel edge-case scenario; returns the shape run.py returns.
%
% One checkpoint per case, keyed by the case name (every name is a valid field name). Each
% case edits a fresh copy of the conventional model and builds a full and a light ecModel
% from it. A model makeEcModel refuses is reported as raised = true with an empty summary;
% the message is not compared. The edit itself sits outside the try block, so a mistake in
% it stops the scenario instead of looking like a refusal.

cases = cellstr(ctx.inputs.cases);
for k = 1:numel(cases)
    results.(cases{k}).full  = build(cases{k}, ctx.inputs, false);
    results.(cases{k}).light = build(cases{k}, ctx.inputs, true);
end

end


function out = build(caseName, inputs, light)
adapter = ModelAdapterManager.getAdapter(inputs.adapter_matlab);
model = loadConventionalGEM('modelAdapter', adapter);
[model, adapter] = edit(caseName, model, adapter);
try
    ecModel = makeEcModel(model, light, adapter);
    out.raised = false;
    out.model = summary(ecModel);
catch
    out.raised = true;
    out.model = struct();
end
end


function [model, adapter] = edit(caseName, model, adapter)
switch caseName
    case 'baseline'
    case 'partial_complex'
        model = changeGrRules(model, {'R3'}, {'G4 and GX'}, true);
    case 'nested_gpr'
        model = changeGrRules(model, {'R2'}, {'G1 and (G2 or G3)'}, true);
    case 'neg_bounds'
        idx = find(strcmp(model.rxns, 'R3'));
        model.lb(idx) = -10;
        model.ub(idx) = -1;
    case 'multi_sink'
        sink.rxns      = {'SINK2'};
        sink.rxnNames  = {'SINK2'};
        sink.equations = {'m1c + m2c <=> '};
        sink.lb        = -1000;
        sink.ub        = 1000;
        model = addRxns(model, sink, 1);
    case 'rxn_id_rev_suffix'
        model.rxns{strcmp(model.rxns, 'R1')} = 'R1_REV';
    case 'rxn_id_exp_suffix'
        model.rxns{strcmp(model.rxns, 'R1')} = 'R1_EXP_1';
    case 'rxn_id_usage_prefix'
        model.rxns{strcmp(model.rxns, 'R5')} = 'usage_prot_P5';
    case 'met_id_prot_prefix'
        model.mets{strcmp(model.mets, 'm2c')} = 'prot_m2';
    case 'enzyme_comp_is_id'
        model.compNames{strcmp(model.comps, 'c')} = 'cytosol';
        model.compNames{strcmp(model.comps, 'e')} = 'extra';
    case 'enzyme_comp_unknown'
        adapter.params.enzyme_comp = 'zz';
    case 'no_genes'
        nGenes = numel(model.genes);
        model.grRules = repmat({''}, numel(model.rxns), 1);
        model.rxnGeneMat = sparse(numel(model.rxns), 0);
        for f = {'genes', 'geneShortNames', 'geneMiriams', 'geneFrom', 'geneComps', 'proteins'}
            if isfield(model, f{1}) && size(model.(f{1}), 1) == nGenes
                model.(f{1}) = model.(f{1})(false(nGenes, 1), :);
            end
        end
    otherwise
        error('unknown case %s', caseName);
end
end


function out = summary(model)
out.gecko_light   = logical(model.ec.geckoLight);
out.n_reactions   = numel(model.rxns);
out.n_metabolites = numel(model.mets);

[rxnIds, order] = sort(model.rxns(:));
reactions = cell(1, numel(rxnIds));
for k = 1:numel(rxnIds)
    i = order(k);
    [lowerValue, lowerKind] = bound(model.lb(i));
    [upperValue, upperKind] = bound(model.ub(i));
    reactions{k} = struct( ...
        'id',                    rxnIds{k}, ...
        'lower_bound',           lowerValue, ...
        'lower_kind',            lowerKind, ...
        'upper_bound',           upperValue, ...
        'upper_kind',            upperKind, ...
        'objective_coefficient', double(model.c(i)));
end
out.reactions = reactions;
out.metabolites = reshape(sort(model.mets(:)), 1, []);
out.stoichiometry = pairs(model.S.', model.rxns, model.mets, {'reaction', 'metabolite'});

ec = model.ec;
out.ec_rxns = reshape(sort(ec.rxns(:)), 1, []);
[enzIds, order] = sort(ec.enzymes(:));
genes = ec.genes(:);
mw = double(ec.mw(:));
enzymes = cell(1, numel(enzIds));
for k = 1:numel(enzIds)
    i = order(k);
    enzymes{k} = struct('enzyme', enzIds{k}, 'gene', genes{i}, 'mw', mw(i));
end
out.ec_enzymes = enzymes;
out.ec_coupling = pairs(ec.rxnEnzMat, ec.rxns, ec.enzymes, {'reaction', 'enzyme'});
end


function [value, kind] = bound(x)
if isinf(x) && x > 0
    value = 0;
    kind = '+inf';
elseif isinf(x) && x < 0
    value = 0;
    kind = '-inf';
else
    value = double(x);
    kind = 'finite';
end
end


function out = pairs(matrix, rowIds, colIds, names)
% Non-zero entries as records (row id, column id, coefficient), sorted on the two keys
% joined with char(1), which reproduces Python's tuple ordering.
[rows, cols, values] = find(matrix);
keys = cell(numel(rows), 1);
for k = 1:numel(rows)
    keys{k} = [rowIds{rows(k)} char(1) colIds{cols(k)}];
end
[~, order] = sort(keys);
out = cell(1, numel(order));
for k = 1:numel(order)
    i = order(k);
    out{k} = struct(names{1}, rowIds{rows(i)}, names{2}, colIds{cols(i)}, ...
        'coefficient', double(values(i)));
end
end
