function results = fseof_toy(ctx)
% MATLAB side of the FSEOF scenario on a hand-written model.

inputs = ctx.inputs;
mets = as_list(inputs.mets);
reactions = as_list(inputs.reactions);
cases = as_list(inputs.cases);
results = struct();
for k = 1:numel(cases)
    c = cases{k};
    model = toy_model(mets, reactions);
    evalc("targets = FSEOF(model, char(inputs.biomass), char(inputs.target), 'iterations', double(c.iterations), 'coefficient', double(c.coefficient), 'corrThreshold', double(c.threshold));");
    out = struct();
    slopes = struct();
    anyTarget = targets.amplify | targets.knockdown | targets.knockout;
    for i = find(anyTarget(:))'
        slopes.(model.rxns{i}) = targets.slope(i);
    end
    out.amplify = sort(model.rxns(targets.amplify))';
    out.knockdown = sort(model.rxns(targets.knockdown))';
    out.knockout = sort(model.rxns(targets.knockout))';
    out.slopes = slopes;
    results.(char(c.name)) = out;
end
end


function model = toy_model(mets, reactions)
n = numel(reactions);
m = numel(mets);
model = struct();
model.id = 'fseof_toy';
model.mets = cell(m, 1);
model.metNames = cell(m, 1);
compartments = cell(m, 1);
for i = 1:m
    model.mets{i} = char(mets{i}.id);
    model.metNames{i} = char(mets{i}.name);
    compartments{i} = char(mets{i}.compartment);
end
[model.comps, ~, model.metComps] = unique(compartments, 'stable');
model.comps = model.comps(:)';
model.compNames = model.comps;
model.metComps = model.metComps(:);
model.rxns = cell(n, 1);
model.rxnNames = cell(n, 1);
model.lb = zeros(n, 1);
model.ub = zeros(n, 1);
model.c = zeros(n, 1);
S = zeros(m, n);
for j = 1:n
    r = reactions{j};
    model.rxns{j} = char(r.id);
    model.rxnNames{j} = char(r.id);
    model.lb(j) = double(r.lb);
    model.ub(j) = double(r.ub);
    model.c(j) = double(r.c);
    for e = as_list(r.stoichiometry)
        S(strcmp(model.mets, char(e{1}.met)), j) = double(e{1}.coefficient);
    end
end
model.S = sparse(S);
model.rev = double(model.lb < 0);
model.b = zeros(m, 1);
model.genes = cell(0, 1);
model.grRules = repmat({''}, n, 1);
model.rxnGeneMat = sparse(n, 0);
end


function out = as_list(value)
% jsondecode gives a struct array when every element has the same fields and a cell otherwise.
if iscell(value)
    out = value(:)';
else
    out = num2cell(value(:)');
end
end
