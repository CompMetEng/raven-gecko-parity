function results = write_fasta_toy(ctx)
% MATLAB side of the writeFasta scenario.

inputs = ctx.inputs;
recs = as_list(inputs.records);
fastaStruct = struct('Header', {}, 'Sequence', {});
for i = 1:numel(recs)
    fastaStruct(i).Header = char(recs{i}.header);
    fastaStruct(i).Sequence = char(recs{i}.sequence);
end
cases = as_list(inputs.cases);
results = struct();
for k = 1:numel(cases)
    c = cases{k};
    file = [tempname() '.fa'];
    writeFasta(file, fastaStruct, double(c.width));
    results.(char(c.name)) = struct('text', fileread(file));
    delete(file);
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
