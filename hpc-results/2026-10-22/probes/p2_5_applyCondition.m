model = pr_common('model'); vals = {'out','in','both','uptake','excrete','blocked','reverse','true','false','none','all','0','1'};
res = struct();
for i = 1:numel(vals)
    f = [tempname '.yml']; fid = fopen(f, 'w'); fprintf(fid, 'prelude:\n  reset_exchanges: %s\n', vals{i}); fclose(fid);
    key = ['v_' vals{i}]; r = struct('yaml', fileread(f));
    try
        evalc('m2 = applyCondition(model, f);');
        ch = find(m2.lb ~= model.lb | m2.ub ~= model.ub);
        r.nChanged = numel(ch); r.ids = m2.rxns(ch); r.lbBefore = model.lb(ch); r.ubBefore = model.ub(ch); r.lbAfter = m2.lb(ch); r.ubAfter = m2.ub(ch);
    catch e, r.error = pr_common('err', e); end
    res.(key) = r; delete(f);
end
pr_common('save', '2.5_applyCondition', res);
