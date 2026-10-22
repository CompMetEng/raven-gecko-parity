function results = uniprot_scores_toy(ctx)
% MATLAB side of the UniProt localisation scenario on a local export.
%
% parseScores with the uniprot predictor is the parser getUniProtScores applies after the download.

inputs = ctx.inputs;
cases = as_list(inputs.cases);
results = struct();
for k = 1:numel(cases)
    c = cases{k};
    args = {'predictor', 'uniprot', 'compartmentMap', compartment_map(inputs, char(c.map))};
    if ~isempty(c.id_column)
        args = [args, {'idColumn', char(inputs.(char(c.id_column)))}]; %#ok<AGROW>
    end
    warning('off', 'RAVEN:warning');
    GSS = parseScores(char(inputs.(char(c.file))), args{:});
    warning('on', 'RAVEN:warning');
    genes = reshape(GSS.genes, 1, []);
    compartments = reshape(GSS.compartments, 1, []);
    byGene = struct();
    for i = 1:numel(genes)
        entry = struct();
        for j = 1:numel(compartments)
            if GSS.scores(i, j) ~= 0
                entry.(compartments{j}) = full(double(GSS.scores(i, j)));
            end
        end
        byGene.(genes{i}) = entry;
    end
    results.(char(c.name)) = struct('genes', {genes}, 'compartments', {compartments}, ...
        'scores', {num2cell(full(double(GSS.scores)), 2)'}, 'by_gene', byGene);
end
end


function map = compartment_map(inputs, name)
if strcmp(name, 'default')
    map = defaultCompartmentMap();
else
    spec = inputs.(name);
    keys = fieldnames(spec);
    values = cellfun(@(f) char(spec.(f)), keys, 'UniformOutput', false);
    map = containers.Map(strrep(keys, '_', ' '), values);
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
