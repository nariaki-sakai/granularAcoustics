function paramSolver = genParamSolverRelax(minPeriodSampling,nIter,stopOpt,stopThreshold,freqSave,f,optBC)
    
    dt0 = 0.01;
    if ~isempty(f)
        dt = floor(min(1/(minPeriodSampling*f),0.1)/dt0)*dt0; %Time step
    else
        dt = 0.01;
    end

    paramSolver.dt = dt;
    paramSolver.nIter = nIter;
    paramSolver.freqSave = freqSave;
    paramSolver.stopOpt = stopOpt;
    paramSolver.stopThreshold = stopThreshold;
    paramSolver.optDisp = 0;
    paramSolver.optOut = 'guide';
    paramSolver.optEflux = 0;
    paramSolver.minIter = nIter;
    paramSolver.maxIter = inf;
    if iscell(optBC)
        paramSolver.optBC = optBC;
    else
        paramSolver.optBC = {optBC;optBC;optBC;optBC};
    end

    
end










