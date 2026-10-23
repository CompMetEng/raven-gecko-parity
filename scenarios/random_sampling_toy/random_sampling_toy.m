function results = random_sampling_toy(ctx)
% MATLAB side of the randomSampling scenario: properties of the sample matrix that do not depend on the generator.

inputs = ctx.inputs;
mets = as_list(inputs.mets);
reactions = as_list(inputs.reactions);
methods = as_list(inputs.methods);
n = double(inputs.n_samples);
results = struct();
for k = 1:numel(methods)
    method = char(methods{k});
    key = method;
    model = toy_model(mets, reactions);
    if strcmp(method, 'random_objective')
        method = 'randomObjective';
    end
    try
        evalc("samples = randomSampling(model, n, 'method', method, 'seed', 1);");
    catch
        results.(key) = struct('error', true);
        continue
    end
    i1 = find(strcmp(model.rxns, 'R1'));
    i2 = find(strcmp(model.rxns, 'R2'));
    mean_flux = mean(samples, 2);
    lb = model.lb;
    ub = model.ub;
    results.(key) = struct( ...
        'n_reactions', size(samples, 1), ...
        'n_samples', size(samples, 2), ...
        'steady_state', max(max(abs(model.S * samples))) < 1e-6, ...
        'steady_state_1e4', max(max(abs(model.S * samples))) < 1e-4, ...
        'within_bounds', all(all(samples >= lb - 1e-6 & samples <= ub + 1e-6)), ...
        'mean_near_centroid', abs(mean_flux(i1) - double(inputs.centroid)) < double(inputs.band), ...
        'parallel_means_close', abs(mean_flux(i1) - mean_flux(i2)) < double(inputs.band));
end
end


function model = toy_model(mets, reactions)
n = numel(reactions);
m = numel(mets);
model = struct();
model.id = 'sampling_toy';
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
