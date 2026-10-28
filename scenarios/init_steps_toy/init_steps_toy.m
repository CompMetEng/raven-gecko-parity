function results = init_steps_toy(ctx)
% MATLAB side of the getINITSteps scenario: the schedule of every series as INITStepDesc objects.

cases = ctx.inputs.cases;
if ~iscell(cases)
    cases = num2cell(cases(:)');
end
results = struct();
for k = 1:numel(cases)
    name = char(cases{k}.series);
    key = char(cases{k}.name);
    try
        steps = getINITSteps('series', name);
    catch
        results.(key) = struct('error', true);
        continue
    end
    records = cell(1, numel(steps));
    for i = 1:numel(steps)
        s = steps{i};
        n = numel(s.MILPParams);
        runs = cell(1, n);
        for j = 1:n
            gap = NaN;
            if iscell(s.AbsMIPGaps)
                gap = double(s.AbsMIPGaps{j});
            elseif ~isempty(s.AbsMIPGaps)
                gap = double(s.AbsMIPGaps(1));
            end
            runs{j} = struct('mip_gap', double(s.MILPParams{j}.MIPGap), 'time_limit', double(s.MILPParams{j}.TimeLimit), 'abs_gap', gap);
        end
        records{i} = struct('how_to_use_prev', s.HowToUsePrevResults, 'ignore_mask', double(s.RxnsToIgnoreMask(:))', ...
            'pos_rev_off', logical(s.PosRevOff), 'allow_met_secr', logical(s.AllowMetSecr), 'runs', {runs});
    end
    results.(key) = struct('steps', {records});
end
end
