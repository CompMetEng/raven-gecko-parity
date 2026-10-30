function results = protein_data_filters_ectestgem(ctx)
% MATLAB side of the protein-data filter scenario.
%
% Returns the shape run.py returns: one entry per case, each {raised, rows}, rows sorted by
% UniProt id with a presence flag and a value per condition. NaN is never emitted (jsonencode
% writes it as null). Parameters are passed by name: a positional call with [] for an
% argument assigns [] literally, which for filterData is false.

adapter = ModelAdapterManager.getAdapter(ctx.inputs.adapter_matlab);

names = sort(fieldnames(ctx.inputs.cases));
results = struct();
for k = 1:numel(names)
    c = ctx.inputs.cases.(names{k});
    file = fullfile(ctx.inputs.data_dir, c.file);
    % A missing table must fail the run: a case that raised because its file is absent would
    % otherwise match the other side's identical failure.
    if ~isfile(file)
        error('protein_data_filters_ectestgem:noFile', 'missing %s', file);
    end
    repl = double(c.repl(:)).';
    maxMissing = double(c.max_missing(:)).';
    try
        protData = loadProtData(repl, ...
            'protDataFile', file, ...
            'filterData',   logical(c.filter_data), ...
            'modelAdapter', adapter, ...
            'minVal',       double(c.min_val), ...
            'maxRSD',       double(c.max_rsd), ...
            'maxMissing',   maxMissing, ...
            'cutLowest',    double(c.cut_lowest), ...
            'addStdevs',    double(c.add_stdevs));
        results.(names{k}) = struct('raised', false, 'rows', {rows_of(protData)});
    catch
        results.(names{k}) = struct('raised', true, 'rows', {{}});
    end
end

end


function out = rows_of(protData)
[ids, order] = sort(protData.uniprotIDs(:));
values = protData.abundances(order, :);

out = cell(1, numel(ids));
for r = 1:numel(ids)
    present = ~isnan(values(r, :));
    v = values(r, :);
    v(~present) = 0;
    out{r} = struct('id', ids{r}, 'present', {num2cell(present)}, 'value', {num2cell(double(v))});
end
end
