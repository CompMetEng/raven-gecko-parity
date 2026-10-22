addpath(getenv('GUROBI_MATLAB'));
res = struct(); res.matlab = version;
for s = {'gurobi','scip','cobra','glpk'}
    r = struct();
    try
        evalc(sprintf('setRavenSolver(''%s'')', s{1}));
        r.setRavenSolver = 'ok'; r.getprefSolver = getpref('RAVEN','solver');
    catch e, r.setRavenSolver = 'error'; r.error = pr_common('err', e); end
    res.(s{1}) = r;
end
res.gurobiOnPath = which('gurobi');
try
    gm.A = sparse([1 1]); gm.obj = [1 1]; gm.rhs = 1; gm.sense = '>'; gm.modelsense = 'min';
    g = gurobi(gm, struct('outputflag', 0)); res.gurobiTinyLP = g.status;
    if isfield(g, 'versioninfo'), res.gurobiVersionInfo = g.versioninfo; end
catch e, res.gurobiTinyLP = pr_common('err', e); end
try, res.gurobiVersion = gurobi_version; catch e, res.gurobiVersion = pr_common('err', e); end
res.grbLicenseFile = getenv('GRB_LICENSE_FILE');
try, res.licenseFileText = strtrim(fileread(getenv('GRB_LICENSE_FILE'))); catch, end
% a tiny MILP on each solver RAVEN offers
for s = {'gurobi','scip','cobra'}
    r = struct();
    try
        evalc(sprintf('setRavenSolver(''%s'')', s{1}));
        p.c = [-1;-1]; p.A = sparse([2 3]); p.b = 7; p.csense = 'L'; p.lb = [0;0]; p.ub = [10;10]; p.vartype = 'II'; p.osense = 1;
        evalc('x = optimizeProb(p, struct(), 1);');
        r.milp = struct('stat', x.stat, 'f', x.obj);
    catch e, r.milpError = pr_common('err', e); end
    res.([s{1} '_milp']) = r;
end
pr_common('save', '3.1_solvers', res);
