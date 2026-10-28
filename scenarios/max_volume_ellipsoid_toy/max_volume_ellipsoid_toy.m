function results = max_volume_ellipsoid_toy(ctx)
% MATLAB side of the sampleMaxVolEllipse scenario.

cases = ctx.inputs.cases;
if ~iscell(cases)
    cases = num2cell(cases(:)');
end
results = struct();
for k = 1:numel(cases)
    c = cases{k};
    A = double(c.A);
    b = double(c.b(:));
    if isempty(c.x0)
        x0 = [];
    else
        x0 = double(c.x0(:));
    end
    [x, E, converged] = sampleMaxVolEllipse(A, b, x0);
    results.(char(c.name)) = struct('centre', x(:)', 'E', {rows(E)}, 'E_Et', {rows(E * E')}, 'converged', logical(converged));
end
end


function out = rows(M)
out = cell(1, size(M, 1));
for i = 1:size(M, 1)
    out{i} = M(i, :);
end
end
