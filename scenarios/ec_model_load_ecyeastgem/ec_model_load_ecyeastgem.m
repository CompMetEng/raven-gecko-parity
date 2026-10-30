function results = ec_model_load_ecyeastgem(ctx)
% MATLAB side of the genome-scale ecModel loading scenario.
%
% loadEcModel reads the ecYeastGEM file of GECKO's full_ecModel tutorial and the loaded
% model is summarised field by field, the way run.py does. See scenario.yml for what is
% compared: lists keep the order of the model, the stoichiometry and the coupling matrix are
% sorted multi-key records joined on char(1), an infinite bound or a NaN is a class and a zero.

adapter = ModelAdapterManager.getAdapter(ctx.inputs.adapter_matlab);
model = loadEcModel(char(ctx.inputs.model_file), adapter);

results.counts = struct( ...
    'reactions', numel(model.rxns), ...
    'metabolites', numel(model.mets), ...
    'genes', numel(model.genes), ...
    'ec_rxns', numel(model.ec.rxns), ...
    'enzymes', numel(model.ec.enzymes));
results.model_id = char(model.id);
results.gecko_light = logical(model.ec.geckoLight);
results.reaction_order = reshape(model.rxns, 1, []);
results.metabolite_order = reshape(model.mets, 1, []);
results.gene_order = reshape(model.genes, 1, []);
results.compartments = reshape(sort(model.comps), 1, []);
results.reactions = reactions(model);
results.gene_rules = gene_rules(model);
results.metabolites = metabolites(model);
results.stoichiometry = pairs(model.S.', model.rxns, model.mets, {'reaction', 'metabolite'}, true);
results.ec = ec_data(model.ec);

end


function out = reactions(model)
n = numel(model.rxns);
out = cell(1, n);
for i = 1:n
    [lowerValue, lowerKind] = bound(model.lb(i));
    [upperValue, upperKind] = bound(model.ub(i));
    out{i} = struct( ...
        'id',                    model.rxns{i}, ...
        'name',                  char(model.rxnNames{i}), ...
        'lower_bound',           lowerValue, ...
        'lower_kind',            lowerKind, ...
        'upper_bound',           upperValue, ...
        'upper_kind',            upperKind, ...
        'objective_coefficient', double(model.c(i)), ...
        'subsystem',             {names(subsystem_of(model, i))}, ...
        'ec_codes',              ec_codes_of(model, i), ...
        'genes',                 {names(model.genes(model.rxnGeneMat(i, :) ~= 0))});
end
end


function out = gene_rules(model)
n = numel(model.rxns);
out = cell(1, n);
for i = 1:n
    rule = '';
    if isfield(model, 'grRules') && numel(model.grRules) >= i
        rule = char(model.grRules{i});
    end
    out{i} = struct('id', model.rxns{i}, 'rule', rule);
end
end


function out = metabolites(model)
n = numel(model.mets);
out = cell(1, n);
for i = 1:n
    charge = NaN;
    if isfield(model, 'metCharges') && numel(model.metCharges) >= i
        charge = double(model.metCharges(i));
    end
    formula = '';
    if isfield(model, 'metFormulas') && numel(model.metFormulas) >= i
        formula = char(model.metFormulas{i});
    end
    out{i} = struct( ...
        'id',          model.mets{i}, ...
        'name',        char(model.metNames{i}), ...
        'compartment', model.comps{model.metComps(i)}, ...
        'formula',     formula, ...
        'charge_set',  ~isnan(charge), ...
        'charge',      zero_for_nan(charge));
end
end


function value = zero_for_nan(x)
value = x;
if isnan(x)
    value = 0;
end
end


function out = names(values)
% A sorted row of the non-empty strings in values, which may be a char, a cell or empty.
if isempty(values)
    out = {};
    return
end
if ischar(values) || isstring(values)
    values = cellstr(values);
end
values = values(:).';
values = values(~cellfun(@isempty, values));
out = sort(values);
end


function value = subsystem_of(model, i)
value = {};
if isfield(model, 'subSystems') && numel(model.subSystems) >= i
    value = model.subSystems{i};
end
end


function value = ec_codes_of(model, i)
value = '';
if isfield(model, 'eccodes') && numel(model.eccodes) >= i && ~isempty(model.eccodes{i})
    value = char(model.eccodes{i});
end
end


function out = ec_data(ec)
out.rxns     = reshape(ec.rxns(:), 1, []);
out.genes    = reshape(ec.genes(:), 1, []);
out.enzymes  = reshape(ec.enzymes(:), 1, []);
[mw, mwKind] = numbers(ec.mw);
out.mw       = mw;
out.mw_kind  = mwKind;
out.sequence = reshape(blanks_for_empty(ec.sequence), 1, []);
[concs, concsKind] = numbers(ec.concs);
out.concs      = concs;
out.concs_kind = concsKind;
out.eccodes  = reshape(blanks_for_empty(ec.eccodes), 1, []);
out.kcat     = num2cell(reshape(double(ec.kcat(:)), 1, []));
out.source   = reshape(blanks_for_empty(ec.source), 1, []);
out.notes    = reshape(blanks_for_empty(ec.notes), 1, []);
out.coupling = pairs(ec.rxnEnzMat, ec.rxns, ec.enzymes, {'reaction', 'enzyme'}, true);
end


function [values, kinds] = numbers(x)
% A number vector as a cell of values and a cell of kinds: NaN and infinity are classes.
x = double(x(:)).';
values = cell(1, numel(x));
kinds = cell(1, numel(x));
for k = 1:numel(x)
    if isnan(x(k))
        values{k} = 0;
        kinds{k} = 'nan';
    else
        [values{k}, kinds{k}] = bound(x(k));
    end
end
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


function out = pairs(matrix, rowIds, colIds, names, withCoefficient)
[rows, cols, values] = find(matrix);
keys = cell(numel(rows), 1);
for k = 1:numel(rows)
    keys{k} = [rowIds{rows(k)} char(1) colIds{cols(k)}];
end
[~, order] = sort(keys);

out = cell(1, numel(order));
for k = 1:numel(order)
    i = order(k);
    record = struct(names{1}, rowIds{rows(i)}, names{2}, colIds{cols(i)});
    if withCoefficient
        record.coefficient = double(values(i));
    end
    out{k} = record;
end
end


function out = blanks_for_empty(values)
out = values(:).';
for k = 1:numel(out)
    if isempty(out{k})
        out{k} = '';
    end
end
end
