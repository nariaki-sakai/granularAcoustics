function paramSolver = genParamSolver(periodSampling,nIterFactor,freqSave,f,paramLattice,optBC)

    %paramSolver = genParamSolver(minPeriodSampling,nIterFactor,freqSave,f,paramLattice,optBC)
    %   generates the structure class variable paramSolver that will
    %   contains all information one would need for numerically solving the
    %   equations of motion
    %
    %   INPUT:
    %       periodSampling: number of sampling during one period of time.
    %           Determines the time step dt according to the excitation
    %           frequency f
    %       nIterFactor: variable that determines the duration of the
    %           simulation. See in the code how it is defined
    %       freqSave: frequency at chich we save the state of the lattice.
    %           freqSave = 1 corresponds to saving every time step
    %       f: central frequency of the source signal. If it does not exist
    %           (for instance with a Dirac), simply put an empty variable
    %           []
    %       paramLattice: the (super useful by its compactness) structure
    %       that contains all information you would need concerning the
    %       lattice (whatever I put inside)
    %       optBC: boundary conditions. needs to be either a string, or a
    %       cell of four strings, where each element is the boundary
    %       condition of one wall (sorted following left wall, bottom, top
    %       and right). Elements should be either 'fixed', 'free' or 'open'
    %       (the latter is for leakage)
    
    
    strain = paramLattice.strain;
    c0 = 2*pi*sqrt(strain/(1+strain));
    
    LAbsLayer = paramLattice.LAbsLayer;
    LtotLayer = paramLattice.LtotLayer;
    L = paramLattice.L;
    Lbuffer = LtotLayer - LAbsLayer;
    
    if ~isempty(f)
        dt = floor(min(1/(periodSampling*f),0.1)/0.01)*0.01; %Time step
    else
        dt = 0.01;
    end
    
    nIter = ceil(nIterFactor*(2*Lbuffer+L(2))/c0/dt/1000)*1000;

    paramSolver.dt = dt;
    paramSolver.nIter = nIter;
    paramSolver.freqSave = freqSave;
    paramSolver.optDisp = 1;
    paramSolver.optOut = 'guide';
    paramSolver.optEflux = 0;
    
    if ischar(optBC)
        paramSolver.optBC = {optBC;optBC;optBC;optBC};
    else
        paramSolver.optBC = optBC;
    end
    
end










