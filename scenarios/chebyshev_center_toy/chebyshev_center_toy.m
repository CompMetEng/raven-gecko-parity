function results = chebyshev_center_toy(ctx)
% MATLAB side of the sampleChebyshevCenter scenario.

cases = ctx.inputs.cases;
if ~iscell(cases)
    cases = num2cell(cases(:)');
end
results = struct();
for k = 1:numel(cases)
    c = cases{k};
    A = double(c.A);
    b = double(c.b(:));
    x0 = sampleChebyshevCenter(A, b);
    radius = min((b - A * x0) ./ sqrt(sum(A.^2, 2)));
    record = struct('radius', radius);
    if c.unique_centre
        record.centre = x0(:)';
    end
    results.(char(c.name)) = record;
end
end
