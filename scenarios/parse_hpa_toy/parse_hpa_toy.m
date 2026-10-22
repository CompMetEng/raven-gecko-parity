function results = parse_hpa_toy(ctx)
% MATLAB side of the parseHPA and parseHPArna scenario: records, independent of the matrix layout.

inputs = ctx.inputs;
results = struct();
protein = as_list(inputs.protein);
for k = 1:numel(protein)
    c = protein{k};
    try
        data = parseHPA(char(c.file));
    catch
        results.(char(c.name)) = struct('error', true);
        continue
    end
    [gi, ti, ~] = find(data.gene2Level);
    records = cell(1, numel(gi));
    for q = 1:numel(gi)
        records{q} = sprintf('%s|%s|%s|%s|%s', data.genes{gi(q)}, data.tissues{ti(q)}, data.celltypes{ti(q)}, ...
            data.levels{data.gene2Level(gi(q), ti(q))}, data.reliabilities{data.gene2Reliability(gi(q), ti(q))});
    end
    results.(char(c.name)) = struct('genes', {sort(data.genes(:))'}, 'records', {sort(records)});
end
rna = as_list(inputs.rna);
for k = 1:numel(rna)
    c = rna{k};
    try
        data = parseHPArna(char(c.file));
    catch
        results.(char(c.name)) = struct('error', true);
        continue
    end
    [gi, ti] = find(~isnan(data.levels));
    records = cell(1, numel(gi));
    for q = 1:numel(gi)
        records{q} = sprintf('%s|%s|%s|%s', data.genes{gi(q)}, data.geneNames{gi(q)}, data.tissues{ti(q)}, ...
            sprintf('%.10g', data.levels(gi(q), ti(q))));
    end
    present = any(~isnan(data.levels), 2);
    results.(char(c.name)) = struct('genes', {sort(data.genes(present))'}, ...
        'tissues', {sort(data.tissues(any(~isnan(data.levels), 1)))'}, 'records', {sort(records)});
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
