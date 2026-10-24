function results = scale_biomass_toy(ctx)
% MATLAB side of the biomass rescaling and fraction scenario on hand-written models.
%
% Each case applies one operation to a fresh model: rescale a pseudoreaction by a factor, scale a component
% to a mass fraction (optionally balancing a second component), or compute the fractions. A case records the
% coefficients of every reaction, and the fractions where the operation produces them; a refused operation is
% recorded as an error.

inputs = ctx.inputs;
cases = as_list(inputs.cases);
results = struct();
for k = 1:numel(cases)
    results.(char(cases{k}.name)) = run_case(inputs, cases{k});
end
end


function out = run_case(inputs, c)
model = build_model(inputs, c);
cfg = build_config(inputs, c);
state = warning('off', 'all');
restoreWarnings = onCleanup(@() warning(state));
try
    switch char(c.op)
        case 'rescale'
            model = scaleBiomassPseudoreaction(model, cfg, char(c.component), double(c.factor));
            out = struct('coefficients', {coefficients(model)});
        case 'scale'
            args = {};
            if ~isempty(c.balance_out)
                args = [args, {'balanceOut', char(c.balance_out)}];
            end
            if ~isempty(c.target_total)
                args = [args, {'targetTotal', double(c.target_total)}];
            end
            model = scaleBiomassFraction(model, cfg, char(c.component), double(c.new_value), args{:});
            out = struct('coefficients', {coefficients(model)}, 'fractions', getBiomassFractions(model, cfg));
        otherwise
            out = struct('fractions', getBiomassFractions(model, cfg));
    end
catch
    out = struct('error', true);
end
end


function model = build_model(inputs, c)
mets = as_list(inputs.models.(char(c.model)));
reactions = as_list(inputs.reactions);
formulaOverrides = as_list(c.formula_overrides);
chargeOverrides = as_list(c.charge_overrides);
m = numel(mets);
n = numel(reactions);
model = struct();
model.id = 'biomass_toy';
model.mets = cell(m, 1);
model.metNames = cell(m, 1);
model.metFormulas = cell(m, 1);
model.metCharges = zeros(m, 1);
for i = 1:m
    model.mets{i} = char(mets{i}.id);
    model.metNames{i} = char(mets{i}.name);
    model.metFormulas{i} = char(mets{i}.formula);
    model.metCharges(i) = double(mets{i}.charge);
end
for i = 1:numel(formulaOverrides)
    model.metFormulas{strcmp(model.mets, char(formulaOverrides{i}.met))} = char(formulaOverrides{i}.formula);
end
for i = 1:numel(chargeOverrides)
    value = chargeOverrides{i}.charge;
    if isempty(value)
        value = NaN;
    end
    model.metCharges(strcmp(model.mets, char(chargeOverrides{i}.met))) = double(value);
end
if c.omit_charges
    model = rmfield(model, 'metCharges');
end
model.comps = {'c'};
model.compNames = {'c'};
model.metComps = ones(m, 1);
model.rxns = cell(n, 1);
model.rxnNames = cell(n, 1);
model.lb = zeros(n, 1);
model.ub = zeros(n, 1);
S = zeros(m, n);
for j = 1:n
    r = reactions{j};
    model.rxns{j} = char(r.id);
    model.rxnNames{j} = char(r.name);
    model.lb(j) = double(r.lb);
    model.ub(j) = double(r.ub);
    for e = as_list(r.stoichiometry)
        S(strcmp(model.mets, char(e{1}.met)), j) = double(e{1}.coefficient);
    end
    if c.proton_in_protein && strcmp(model.rxns{j}, 'P1')
        S(strcmp(model.mets, 'h_c'), j) = -0.1;
    end
end
model.S = sparse(S);
model.rev = double(model.lb < 0);
model.c = zeros(n, 1);
model.b = zeros(m, 1);
model.genes = cell(0, 1);
model.grRules = repmat({''}, n, 1);
model.rxnGeneMat = sparse(n, 0);
end


function cfg = build_config(inputs, c)
components = as_list(inputs.configs.(char(c.config)));
cfg = struct();
cfg.biomass_rxn = 'P1';
cfg.proton_met = char(c.proton_met);
cfg.components = cell(1, numel(components));
for i = 1:numel(components)
    cfg.components{i} = struct('name', char(components{i}.name), ...
        'pseudoreaction_name', char(components{i}.pseudoreaction_name), ...
        'mass_strategy', char(components{i}.mass_strategy));
end
end


function entries = coefficients(model)
entries = {};
keys = {};
for j = 1:numel(model.rxns)
    rows = find(model.S(:, j));
    for q = 1:numel(rows)
        entries{end + 1} = struct('reaction', model.rxns{j}, 'met', model.mets{rows(q)}, ...
            'coefficient', full(double(model.S(rows(q), j)))); %#ok<AGROW>
        keys{end + 1} = [model.rxns{j} '|' model.mets{rows(q)}]; %#ok<AGROW>
    end
end
[~, order] = sort(keys);
entries = entries(order);
end


function out = as_list(value)
% jsondecode gives a struct array when every element has the same fields and a cell otherwise.
if iscell(value)
    out = value(:)';
elseif isempty(value)
    out = {};
else
    out = num2cell(value(:)');
end
end
