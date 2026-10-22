function results = set_gam_toy(ctx)
% MATLAB side of the setGAM scenario on a hand-written model.
%
% Each case sets the cofactors of a biomass reaction to one or more values on a fresh copy of the model and
% records every coefficient of the biomass reaction, whether its magnitude equals the value exactly, and the
% bounds of the maintenance reaction.

inputs = ctx.inputs;
mets = as_list(inputs.mets);
cases = as_list(inputs.cases);
results = struct();
for k = 1:numel(cases)
    c = cases{k};
    runs = cell(1, numel(c.values));
    values = double(c.values(:))';
    for v = 1:numel(values)
        runs{v} = run_one(mets, c, values(v));
    end
    results.(char(c.name)) = struct('runs', {runs});
end
end


function out = run_one(mets, c, value)
model = toy_model(mets, as_list(c.reactions));
cofactors = as_cellstr(c.cofactors);
args = {};
if ~isempty(c.ngam_rxn)
    args = [args, {'ngamRxn', char(c.ngam_rxn)}]; %#ok<AGROW>
end
if ~isempty(c.ngam_value)
    args = [args, {'ngamValue', double(c.ngam_value)}]; %#ok<AGROW>
end
state = warning('off', 'all');
restoreWarnings = onCleanup(@() warning(state));
try
    model = setGAM(model, value, char(c.biomass_rxn), cofactors, args{:});
catch
    out = struct('value', value, 'error', true);
    return
end
bio = find(strcmp(model.rxns, char(c.biomass_rxn)));
rows = find(model.S(:, bio));
entries = cell(1, numel(rows));
for q = 1:numel(rows)
    coef = full(double(model.S(rows(q), bio)));
    entries{q} = struct('met', model.mets{rows(q)}, 'coefficient', coef, 'abs_equals_value', abs(coef) == abs(value));
end
[~, order] = sort(cellfun(@(e) e.met, entries, 'UniformOutput', false));
entries = entries(order);
atpm = strcmp(model.rxns, 'ATPM');
out = struct('value', value, 'coefficients', {entries}, 'ngam_lb', model.lb(atpm), 'ngam_ub', model.ub(atpm));
end


function model = toy_model(mets, reactions)
n = numel(reactions);
m = numel(mets);
model = struct();
model.id = 'gam_toy';
model.mets = cell(m, 1);
model.metNames = cell(m, 1);
compartments = {};
for i = 1:m
    model.mets{i} = char(mets{i}.id);
    model.metNames{i} = char(mets{i}.name);
    compartments{i} = char(mets{i}.compartment); %#ok<AGROW>
end
[model.comps, ~, model.metComps] = unique(compartments(:), 'stable');
model.comps = model.comps(:)';
model.compNames = model.comps;
model.metComps = model.metComps(:);
model.rxns = cell(n, 1);
model.rxnNames = cell(n, 1);
model.lb = zeros(n, 1);
model.ub = zeros(n, 1);
S = zeros(m, n);
for j = 1:n
    r = reactions{j};
    model.rxns{j} = char(r.id);
    model.rxnNames{j} = char(r.id);
    model.lb(j) = double(r.lb);
    model.ub(j) = double(r.ub);
    for e = as_list(r.stoichiometry)
        S(strcmp(model.mets, char(e{1}.met)), j) = double(e{1}.coefficient);
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
