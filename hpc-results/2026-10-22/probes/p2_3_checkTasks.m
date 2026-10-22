tc = RavenTestCase_helper();
model = tc;
hdr = sprintf('\tID\tDESCRIPTION\tSHOULD FAIL\tIN\tIN LB\tIN UB\tOUT\tOUT LB\tOUT UB\tEQU\tEQU LB\tEQU UB\tCHANGED RXN\tCHANGED LB\tCHANGED UB\tPRINT FLUX\tCOMMENTS\tREFERENCES\tEXPLANATIONS');
cases = struct();
%            ID  DESC SF IN INLB INUB OUT OUTLB OUTUB EQU EQULB EQUUB CRXN CLB CUB
cases.A = {'a[s]','','','e[s]','1','','','','','NOSUCHRXN','0','0'};
cases.B = {'a[s]','2','','a[s];e[s]','3','','','','','','',''};
cases.C = {'a[s]; a[s]','','','e[s]','1','','','','','','',''};
cases.D = {'a[s]','','','e[s]','1','','','','','R2','5','1'};
cases.E = {'a[s]','','','e[s]','1','','','','','R2; R3','-5; 10',''};
cases.F = {'a[s]','','','e[s]','1','','a[c] => b[c]','7','3','','',''};
res = struct();
for k = fieldnames(cases)'
    n = k{1}; c = cases.(n);
    f = [tempname '.txt']; fid = fopen(f, 'w'); fprintf(fid, '%s\n', hdr);
    fprintf(fid, '\t%s\t%s\t0\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t0\t\t\t\n', n, n, c{:});
    fclose(fid);
    r = struct('file', fileread(f));
    try
        evalc('tr = checkTasks(model, f, true, false, false);');
        r.ok = tr.ok; if isfield(tr, 'ok'), r.okAll = all(tr.ok); end
    catch e, r.error = pr_common('err', e); end
    res.(n) = r; delete(f);
end
pr_common('save', '2.3_checkTasks', res);
function m = RavenTestCase_helper()
m = struct();
m.id = 'testModel'; m.rxns = {}; m.S = []; m.rev = [];
m.mets = {'as';'ac';'bc';'cc';'dc';'ec';'es';'fc'};
m.metNames = {'a';'a';'b';'c';'d';'e';'e';'f'};
m.comps = {'s';'c'}; m.compNames = m.comps;
m.metComps = [1;2;2;2;2;2;1;2];
m.genes = {'G1';'G2';'G3';'G4';'G5';'G6';'G7';'G8';'G9';'G10'};
m.grRules = {}; m.rxnGeneMat = [];
r = struct();
r.rxns = {'R1';'R2';'R3';'R4';'R5';'R6';'R7';'R8';'R9';'R10'};
r.equations = {'=> a[s]';'a[s] <=> a[c]';'a[c] <=> b[c] + c[c]';'a[c] <=> 2 d[c]';'b[c] + c[c] => e[c]';'2 d[c] => e[c]';'e[c] => e[s]';'e[s] =>';'a[c] <=> f[c]';'f[c] <=> e[c]'};
r.grRules = {'';'';'G3';'G4';'G5';'G6';'G7';'';'G9';'G10'};
evalc('m = addRxns(m, r, 3);');
m.c = [0;0;0;0;0;0;0;1;0;0]; m.ub = repmat(1000,10,1);
m.lb = [0;-1000;-1000;-1000;0;0;0;0;-1000;-1000];
m.rxnNames = m.rxns; m.b = repmat(0,8,1); m = closeModel(m);
end
