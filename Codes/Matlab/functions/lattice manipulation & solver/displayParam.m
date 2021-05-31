function displayParam(paramLattice,f,paramSolver)

    %displayParam(paramLattice,f,paramSolver)
    %   Displays some useful information about the lattice and the solver
    %   parameters. Useful to remind to the user what he is doing

    LtotLayer = paramLattice.LtotLayer;
    LAbsLayer = paramLattice.LAbsLayer;
    Lbuffer = LtotLayer - LAbsLayer;
    L = paramLattice.L;
    Ltot = [L(1) L(2)+2*LtotLayer];
        
    dt = paramSolver.dt;
    nIter = paramSolver.nIter;
    
    strain = paramLattice.strain;
    c0 = 2*pi*sqrt(strain/(1+strain));
    
    %DISPLAY PARAM
    fprintf('sample size: %d-%d\n',L(1),L(2));
    fprintf('buffer size: %d\n',Lbuffer);
    disp(['Total nb nodes = ' num2str(prod(Ltot))])

    disp(['wave velocity (a/T0) = ' num2str(c0)]);
    try
        lambda0 = c0/f;
        nc = f*L(1)/c0;
        disp(['wave length (a) = ' num2str(lambda0)]);
        disp(['max mode index = ' num2str(2*nc)]);
    catch
    end
    
    disp(['nb iterations = ' num2str(nIter)]);
    disp(['dt = ' num2str(dt)]);

end











