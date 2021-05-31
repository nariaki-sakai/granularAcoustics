function dataOut = verletSolverSH(Ui,dUi,U0,paramLattice,paramSolver,paramSource)


    %dataOut = verletSolverSH(Ui,dUi,U0,paramLattice,paramSolver,paramSource)
    %   solve numerically the equation of motions of all masses for SH out-of-plane
    %   motions. The solver uses the velocty Verlet algorithm
    %
    %   INPUT:
    %       Ui: position of all mass at t=0. Has size (Ly,Lx)
    %       dUi: velocity at t=0. Has size(Lx,Ly) as well
    %       U0: in-plane displacement field at rest
    %       paramLattice: structure that contrains all parameters concerning the
    %           lattice
    %       paramSolver:  structure that contrains all parameters concerning the
    %           numerical solver
    %       paramLSource: structure that contrains all parameters concerning the
    %           source
    %
    %   OUTPUT:
    %       dataOut: structure that contrains all output
    %           dataOut.U: position in time of the masses. We record only a
    %               portion of the system, according to the variable
    %               paramSolver.optOut. Except if paramSolver.optOut =
    %               'detector', dataOut.U has size (nt,Li,Lj) where nt is
    %               the number of recorded time steps, and Li=Ly and Lj are
    %               the widths and length of the portion of the wave guide
    %               that has been recorded. See genParamSolver
    %           dataOut.dU: velocity in time of masses. Same than
    %               dataOut.U, all masses are not recorded, see above
    %               dataOut.time = time;
    %           dataOut.E: kinetic and potential energy of the system VS time
    %           dataOut.Ulast = position of ALL masses (including the absorbing layer)
    %               at the end of the simulation. Useful if one wants to
    %               continue the simulation with this final state, simply
    %               by putting this state to the input of verletSolverSH
    %           dataOut.dUlast: same than above but with the velocity
    %           dataOut.flux: energy exchange of each mass with its
    %           neighbours. has size (nt, Li, Lj, 8), where the last
    %           dimension corresponds to the neighbour, sorted from
    %           - bottom left, middle left, top left (1 2 3)
    %           - bottom middle, top middle (4 5)
    %           - bottom right, middle right, top right
    %           It means in the case there is four neighbours on left right
    %           top bottom, only the elements (2, 4, 5, 7) are non zero
    %             
    %
    
    %%%%%%%%%%%
    %PARAMETERS
    %%%%%%%%%%%
    mMat = paramLattice.mMat;
    if var(mMat(:)) == 0
        mMat = mean(mMat,'all');
    else
    end
    
    nNN = paramLattice.nNN;
    
    dt = paramSolver.dt;
    nIter = paramSolver.nIter;
    freqSave = paramSolver.freqSave;
    nt = floor(nIter/freqSave);
    optDisp = paramSolver.optDisp;    
    
    
    %%%%%%%
    %OUTPUT
    %%%%%%%
    UvsTime = formatVarOut2(Ui,paramLattice,paramSolver.optOut);
    UvsTime = zeros([nt size(UvsTime)]);
    sizeU = size(UvsTime);
    
    dUvsTime = zeros(sizeU);
    if paramSolver.optEflux
        Eflux = nan([sizeU nNN]);
    end
    
    EvsTime = zeros(nt,2);
    time = (1:nt)'*freqSave*dt;
    
    %%%%%%%%%%%
    %INITIALIZE
    %%%%%%%%%%%
    Uz = Ui;
    dUz = dUi;
    Ux = U0(:,:,1);
    Uy = U0(:,:,2);
    Unx = neighbourVal(Ux);
    Uny = neighbourVal(Uy);
    [distNodeX, distNodeY] = distanceNodeXY(Ux,Uy,Unx,Uny);

    [Uz,dUz] = fixBoundaryConditions(Uz,dUz,paramSolver);
    if strcmp(paramSource.forcing,'position')
        [Uz,dUz] = fixSources(Uz,dUz,paramSource,paramSolver,0);
    end
    
    Unz  = neighbourVal(Uz);
    dUnz = neighbourVal(dUz);
    [FzC, ~ ,~] = forceConserv(Uz,Unz,distNodeX,distNodeY,paramLattice);
    FzD = forceDissip(dUz,dUnz,paramLattice);
    d2Uz1 = (FzC + FzD)./mMat;
    d2Uz1(mMat == 0) = 0;
    
    if isnan(FzC)
        error('indt = 0: forceConserv computes nan or inf');
    elseif isnan(FzD)
        error('indt = 0: forceDissip computes nan or inf');
    end
    
    %debugging
    if 0 == 1
        %%
        figure(13);
        subplot(1,3,1);
        cla;
        sourceAmp = paramSource.sourceAmp;
        imagesc(Uz,sourceAmp*[-1 1]);
        axis equal;
        xlim([0 L(2)]+0.5);
        ylim([0 L(1)]+0.5);
                
        subplot(1,3,2);
        cla;
        imagesc(FzC);
        axis equal;
        xlim([0 L(2)]+0.5);
        ylim([0 L(1)]+0.5);
        
        subplot(1,3,3);
        cla;
        imagesc(mMat);
        axis equal;
        xlim([0 L(2)]+0.5);
        ylim([0 L(1)]+0.5);
        
        formatfig(13,[600 200]);
    end
    
    %%
    for indt = 1:nIter

        if optDisp
            dispIndForLoop(indt,nIter,10^(round(log10(nIter/10))));
        end
        %% VERLET
        %POSITION t+1
        Uz = Uz + dUz*dt + d2Uz1*dt^2/2;

        %VELOCITY t+1 PASS 0
        tmpdUz = dUz + d2Uz1*dt;
        
        [Uz,tmpdUz] = fixBoundaryConditions(Uz,tmpdUz,paramSolver);
        if strcmp(paramSource.forcing,'position')
            [Uz,tmpdUz] = fixSources(Uz,tmpdUz,paramSource,paramSolver,indt);
        end
                
        %ACCELERATION AT T+1 PASS 1
        Unz  = neighbourVal(Uz);
        dUnz = neighbourVal(tmpdUz);
        
        [FzC, distNode, FzCAll] = forceConserv(Uz,Unz,distNodeX,distNodeY,paramLattice);
        FzD = forceDissip(tmpdUz,dUnz,paramLattice);
        d2Uz2 = (FzC + FzD)./mMat;
        d2Uz2(mMat == 0) = 0;
        if isnan(FzC)
            error('indt = %d pass 1: forceConserv computes nan or inf',indt);
        elseif isnan(FzD)
            error('indt = %d pass 1: forceDissip computes nan or inf',indt);
        end

        if strcmp(paramSource.forcing,'force')
            d2Uz2 = addForcing(d2Uz2,paramSource,indt);
        end

        %VELOCITY PASS 1
        tmpdUz = dUz + (d2Uz1+d2Uz2)*dt/2;
        
        [Uz,tmpdUz] = fixBoundaryConditions(Uz,tmpdUz,paramSolver);
        if strcmp(paramSource.forcing,'position')
            [Uz,tmpdUz] = fixSources(Uz,tmpdUz,paramSource,paramSolver,indt);
        end
        
        %ACCELERATION AT T+1 PASS 2
        dUnz = neighbourVal(tmpdUz);
        [FzD, FzFAll, FzV] = forceDissip(tmpdUz,dUnz,paramLattice);
        d2Uz2 = (FzC + FzD)./mMat;
        d2Uz2(mMat == 0) = 0;
        if isnan(FzD)
            error('indt = %d pass 2: forceDissip computes nan or inf',indt);
        end
        
        if strcmp(paramSource.forcing,'force')
            d2Uz2 = addForcing(d2Uz2,paramSource,indt);
        end
        
        %VELOCITY PASS 2
        dUz = dUz + (d2Uz1+d2Uz2)*dt/2;
        
        [Uz,dUz] = fixBoundaryConditions(Uz,dUz,paramSolver);
        if strcmp(paramSource.forcing,'position')
            [Uz,dUz] = fixSources(Uz,dUz,paramSource,paramSolver,indt);
        end
        
        %%NEXT
        d2Uz1 = d2Uz2;

        %debugging
        if 0 == 1
            %%
            figure(13);
            cla;
            sourceAmp = paramSource.sourceAmp;
            imagesc(Uz,sourceAmp*[-1 1]);
            
        end
        
        %% Save
        if mod(indt,freqSave) == 0
            
            if paramSolver.optEflux
                flux = energyFlux(dUz,dUnz,FzCAll,FzFAll,FzV,paramLattice);
                Eflux(floor(indt/freqSave),:,:,:)  = formatVarOut3(flux,paramLattice,paramSolver.optOut);
            end
            [Ec, Ep] = energy(distNode,dUz,paramLattice);
            
            UvsTime(floor(indt/freqSave),:,:)  = formatVarOut2(Uz,paramLattice,paramSolver.optOut);
            dUvsTime(floor(indt/freqSave),:,:) = formatVarOut2(dUz,paramLattice,paramSolver.optOut);
            EvsTime(floor(indt/freqSave),1) = sum(Ec(:));
            EvsTime(floor(indt/freqSave),2) = sum(Ep(:));

        end
        
    end
    
    UvsTime = squeeze(UvsTime);
    dUvsTime = squeeze(dUvsTime);
    
    dataOut.U = UvsTime;
    dataOut.dU = dUvsTime;
    dataOut.time = time;
    dataOut.E = EvsTime;
    dataOut.Ulast = Uz;
    dataOut.dUlast = dUz;
    
    if paramSolver.optEflux
        Eflux = squeeze(Eflux);
        dataOut.flux = Eflux;
    end
    
    
    
end

function Un = neighbourVal(U)

    %%
    vecLattice = [-1 -1;
           0 -1;
           1 -1;
          -1  0;
           1  0;
          -1  1;
           0  1;
           1  1];

    %%% X
    Un = zeros([size(U) 8]);
    for ind = 1:8
        Un(:,:,ind) = circshift(U,-vecLattice(ind,:));
    end
    
    %Remove periodic boundaries
    L = size(Un);
    Un(:   ,1,1) = 0;
    Un(:   ,1,2) = 0;
    Un(:   ,1,3) = 0;
    
    Un(:,L(2),6) = 0;
    Un(:,L(2),7) = 0;
    Un(:,L(2),8) = 0;
    
    Un(1,:,1) = 0;
    Un(1,:,4) = 0;
    Un(1,:,6) = 0;
    
    Un(L(1),:,3) = 0;
    Un(L(1),:,5) = 0;
    Un(L(1),:,8) = 0;
    

    
    

end

function distNodeZ = distanceNodeZ(Uz,Unz)

    distNodeZ = zeros(size(Unz));
    
    vecLattice = [-1 -1;
               0 -1;
               1 -1;
              -1  0;
               1  0;
              -1  1;
               0  1;
               1  1];
           
    for indNN = 1:4
        distNodeZ(:,:,indNN) = Unz(:,:,indNN) - Uz;
    end
    
    distNodeZ(:,1,1) = 0;
    distNodeZ(:,1,2) = 0;
    distNodeZ(:,1,3) = 0;
    distNodeZ(1,:,4) = 0;
    
    for indNN = 5:8
        distNodeZ(:,:,indNN) = -circshift(distNodeZ(:,:,9-indNN),vecLattice(9-indNN,:));
    end
        
end

function [distNodeX, distNodeY] = distanceNodeXY(Ux,Uy,Unx,Uny)

    distNodeX = zeros(size(Unx));
    distNodeY = zeros(size(Uny));
    
    vecLattice = [-1 -1;
               0 -1;
               1 -1;
              -1  0;
               1  0;
              -1  1;
               0  1;
               1  1];
    for indNN = 1:4
        distNodeX(:,:,indNN) = Unx(:,:,indNN) + vecLattice(indNN,2) - Ux;
        distNodeY(:,:,indNN) = Uny(:,:,indNN) + vecLattice(indNN,1) - Uy;
    end
    for indNN = 5:8
        distNodeX(:,:,indNN) = -circshift(distNodeX(:,:,9-indNN),vecLattice(9-indNN,:));
        distNodeY(:,:,indNN) = -circshift(distNodeY(:,:,9-indNN),vecLattice(9-indNN,:));
    end
        
end

function [FzC, distNode, FzCAll] = forceConserv(Uz,Unz,distNodeX,distNodeY,paramLattice)
    %%
%     xc = 1/20;
%     xLJ = 2^(1/6);
    
    %%PARAMETERS
    kMat = paramLattice.kMat;
    mMat = paramLattice.mMat;
    a = paramLattice.a;

    %%DIST NODES
    distNodeZ = distanceNodeZ(Uz,Unz);
    distNode = sqrt(distNodeX.^2 + distNodeY.^2 + distNodeZ.^2);
        
    %%Z
    FzCAll = kMat.*(1-a./distNode).*distNodeZ;
    
    
%     if paramLattice.nNN == 4
%         FzCAll(:,:,[1 3 6 8]) == 0;
%     end
%     FzCAll = FzCAll + 
    FzC = sum(FzCAll,3);
    FzC(mMat == 0) = 0;
    hasError = dispNaNInf(FzC);
    if hasError
        warning('Cons. forces has NaN or Inf');
        FzC = nan;
    end
    
    if 0 == 1
        %%
        L = size(Uz);
        sourceAmp = paramSource.sourceAmp;
        
        figure(13);
        subplot(1,3,1);
        cla;
        imagesc(Uz,sourceAmp*[-1 1]);
        axis equal;
        xlim([0 L(2)]+0.5);
        ylim([0 L(2)]+0.5);
        subplot(1,3,2);
        cla;
        tmp = kMat;
        imagesc(tmp(:,:,7));
        axis equal;
        xlim([0 L(2)]+0.5);
        ylim([0 L(2)]+0.5);
        subplot(1,3,3);
        cla;
        imagesc(mMat);
        axis equal;
        xlim([0 L(2)]+0.5);
        ylim([0 L(2)]+0.5);
        formatfig(13,[600 200]);

    end
        
end

function [FzD, FzFAll, FzV] = forceDissip(dUz,dUnz,paramLattice)

    %%PARAMETERS
    viscousMat = paramLattice.viscous;
    frictionMat = paramLattice.friction;
    mMat = paramLattice.mMat;

    %%VISCOUS
    FzV = - viscousMat .* dUz;
    FzV(mMat == 0) = 0;
    
    %%% FRICTION X
    dUz0 = reshape(dUz,[size(dUz) 1]);
    relVelZ = dUnz - dUz0;
    FzFAll = frictionMat.*relVelZ;
    FzF = sum(FzFAll,3);
    FzF(mMat == 0) = 0;

    %%% CHECK AND OUTPUT
    hasError1 = dispNaNInf(FzF);
    hasError2 = dispNaNInf(FzV);
    FzD = FzF + FzV;
    if hasError1*hasError2
        warning('Dissp. forces has NaN or Inf');
        FzD = nan;
    end

end

function [Uz,dUz] = fixSources(Uz,dUz,paramSource,paramSolver,indt)
    
%%
    dt = paramSolver.dt;
    try 
        USources = paramSource.USources;
        sourcePos = paramSource.sourcePos;
        ntSource = size(USources,1);
        nSources = size(USources,2);
    catch
        ntSource = 0;
    end
    
    indt = indt + 1;
    if indt <= ntSource
        for indSource=1:nSources
            indij = sourcePos(indSource,:);
            indi = indij(1);
            indj = indij(2);

            Uz(indi,indj)  =  USources(indt,indSource);
            if indt > 1 && indt < ntSource
                dUz(indi,indj)  =  (USources(indt+1,indSource,1)-USources(indt-1,indSource,1))/(2*dt);
            elseif indt == 1
                dUz(indi,indj)  =  (USources(indt+1,indSource,1)-USources(indt,indSource,1))/dt;
            else
                dUz(indi,indj)  =  (USources(indt,indSource,1)-USources(indt-1,indSource,1))/dt;
            end
        end
    end
    
end

function [Uz,dUz] = fixBoundaryConditions(Uz,dUz,paramSolver)
    
%     try 
%         wallPos = paramSolver.wallPos;
%         wallVel = paramSolver.wallVel;
%     catch
%         optBC = paramSolver.optBC;
%         switch optBC
%             case 'fixed'
%                 wallPos = zeros(4,1);
%                 wallVel = zeros(4,2);
%             case 'free'
%             otherwise
%                 error('paramSolver.optBC has wrong attribute');
%         end
%     end
    
    optBC = paramSolver.optBC;
    
    
    L = size(Uz);
    
    indij = cell(4,2);
    indij{1,1} = 1:L(1);
    indij{1,2} = 1;
    indij{2,1} = 1;
    indij{2,2} = 1:L(2);
    indij{3,1} = L(1);
    indij{3,2} = 1:L(2);
    indij{4,1} = 1:L(1);
    indij{4,2} = L(2);
    
%     for indB = 1:4
%         if ~isnan(wallPos(indB))
%             indi = indij{indB,1};
%             indj = indij{indB,2};
%             Uz(indi,indj) = wallPos(indB);
%         end
%         if ~isnan(wallVel(indB))
%             indi = indij{indB,1};
%             indj = indij{indB,2};
%             dUz(indi,indj) = wallVel(indB);
%         end
%     end

    for indB = 1:4
        tmpBC = optBC{indB};
        switch class(tmpBC)
            case 'char'
                switch tmpBC
                    case 'clamped'
                        wallPos = zeros(4,1);
                        wallVel = zeros(4,2);
                        indi = indij{indB,1};
                        indj = indij{indB,2};
                        Uz(indi,indj) = wallPos(indB);
                        
                        indi = indij{indB,1};
                        indj = indij{indB,2};
                        dUz(indi,indj) = wallVel(indB);
                    case 'free'
                    case 'open'
                    otherwise
                        error('paramSolver.optBC{%s} has wrong attribute',indBC);
                end
            case 'double'
            otherwise
                error('paramSolver.optBC{%s} has data type mismatch',indBC);
        end
    end
    
    
end

function d2Uz = addForcing(d2Uz,paramSource,indt)

    m0 = 1/(2*pi)^2;
    try 
        USources = paramSource.USources;
        sourcePos = paramSource.sourcePos;
        ntSource = size(USources,1);
        nSources = size(USources,2);
    catch
        ntSource = 0;
    end
    if indt <= ntSource
        for indSource=1:nSources
            indij = sourcePos(indSource,:);
            indi = indij(1);
            indj = indij(2);

            d2Uz(indi,indj) = d2Uz(indi,indj) + USources(indt,indSource)/m0;
        end
    end

end

function Eflux = energyFlux(dUz,dUnz,FzCAll,FzFAll,FzV,paramLattice)

    %%
    %%%INPUT
    Ltot = size(dUz);
    try
        nNN = paramLattice.nNN;
    catch
        nNN = 4;
    end

    %%%OUTPUT
    Eflux = zeros(Ltot(1),Ltot(2),nNN);

    %%%FLUX
    if nNN == 4
        
        for ind = 1:4
            indNN = 2*ind;
            if indNN > 5
                indNN = indNN - 1;
            end
            Eflux(:,:,ind) = FzCAll(:,:,indNN).*(dUnz(:,:,indNN,1)+dUz(:,:,1))/2 + ...
                               FzFAll(:,:,indNN).*dUz(:,:,1) + ...
                               FzV(:,:,1).*dUz(:,:,1);

        end
        
    elseif nNN == 8
        for indNN = 1:8
            Eflux(:,:,indNN) = FzCAll(:,:,indNN).*(dUnz(:,:,indNN,1)+dUz(:,:,1))/2 + ...
                                 FzFAll(:,:,indNN).*dUz(:,:,1) + ...
                                 FzV(:,:,1).*dUz(:,:,1);

        end
        
    else
        error('number of neighbours should be 4 or 8');
    end

end

function varOut = formatVarOut2(varIn,paramLattice,optOut)

    L = paramLattice.L;
    Lbuffer = paramLattice.Lbuffer;
    LAbsLayer = paramLattice.LAbsLayer;
    LtotLayer = paramLattice.LtotLayer;
    Ltot = [L(1) 2*LtotLayer+L(2)];
    
    switch optOut
        case 'full'
            varOut  = varIn;
        case 'part'
            if 2*Lbuffer+3*L(2) < Ltot(2) - LAbsLayer
                try
                    switch paramLattice.removeLead
                        case 'left'
                            varOut = varIn(:,1:LtotLayer+2*L(2));
                        case 'right'
                            varOut = varIn(:,LAbsLayer+1-L(2):end);
                    end
                catch
                    varOut = varIn(:,LAbsLayer+1-L(2):LtotLayer+2*L(2)+Lbuffer);
                end
            else
                try
                    switch paramLattice.removeLead
                        case 'left'
                            varOut = varIn(:,1:end-LAbsLayer/2);
                        case 'right'
                            varOut = varIn(:,LAbsLayer/2+1:end);
                    end
                catch
                    varOut = varIn(:,LAbsLayer/2+1:end-LAbsLayer/2);
                end
            end                    
        case 'guide'
            try
                switch paramLattice.removeLead
                    case 'left'
                        varOut  = varIn(:,1:Lbuffer+L(2));
                    case 'right'
                        varOut  = varIn(:,LAbsLayer+1:LtotLayer+L(2));
                end
            catch
               varOut  = varIn(:,LAbsLayer+1:LtotLayer+L(2)+Lbuffer);
            end
        case 'sample'
            try
                switch paramLattice.removeLead
                    case 'left'
                        varOut  = varIn(:,1:L(2));
                    case 'right'
                        varOut  = varIn(:,LtotLayer+1:LtotLayer+L(2));
                end
            catch
                varOut  = varIn(:,LtotLayer+1:LtotLayer+L(2));
            end
        case 'detector'
            detectorPos = paramLattice.detectorPos;
            detectorPos = cell2mat(detectorPos);
            indDetector = sub2ind(Ltot,detectorPos(:,1),detectorPos(:,2));
            indDetector = sort(indDetector,'ascend');
            indDetector = unique(indDetector);
            nDetector = length(indDetector);
            
            varOut = zeros(1,nDetector);
            for indD = 1:nDetector
                [indi,indj] = ind2sub(Ltot,indDetector(indD));
                varOut(1,indD) = varIn(indi,indj);
            end
        otherwise
            error('paramSolver.fullOutput: wrong attribute, needs to be among full, sample or detector');
    end

end

function varOut = formatVarOut3(varIn,paramLattice,optOut)

    L = paramLattice.L;
    Lbuffer = paramLattice.Lbuffer;
    LAbsLayer = paramLattice.LAbsLayer;
    LtotLayer = paramLattice.LtotLayer;
    Ltot = [L(1) 2*LtotLayer+L(2)];
    nNN = paramLattice.nNN;
    
    switch optOut
        case 'full'
            varOut  = varIn;
        case 'part'
            if 2*Lbuffer+3*L(2) < Ltot(2) - LAbsLayer
                varOut = varIn(:,LAbsLayer+1-L(2):LtotLayer+2*L(2)+Lbuffer,:);
            else
                varOut = varIn(:,LAbsLayer/2+1:end-LAbsLayer/2,:);
            end                    
        case 'guide'
            varOut  = varIn(:,LAbsLayer+1:LtotLayer+L(2)+Lbuffer,:);
        case 'sample'
            varOut  = varIn(:,LtotLayer+1:LtotLayer+L(2),:);
        case 'detector'
            detectorPos = paramLattice.detectorPos;
            detectorPos = cell2mat(detectorPos);
            indDetector = sub2ind(Ltot,detectorPos(:,1),detectorPos(:,2));
            indDetector = sort(indDetector,'ascend');
            indDetector = unique(indDetector);
            nDetector = length(indDetector);
            
            varOut = zeros(1,nDetector,nNN);
            for indD = 1:nDetector
                [indi,indj] = ind2sub(Ltot,indDetector(indD));
                varOut(1,indD,:) = varIn(indi,indj,:);
            end
        otherwise
            error('paramSolver.fullOutput: wrong attribute, needs to be among full, sample or detector');
    end

end

function varOut = formatRemoveLead(varIn,paramLattice,optOut)

    

end

% function varOut = formatVarOut4(varIn,paramLattice,optOut)
% 
%     L = paramLattice.L;
%     Lbuffer = paramLattice.Lbuffer;
%     LAbsLayer = paramLattice.LAbsLayer;
%     LtotLayer = paramLattice.LtotLayer;
%     Ltot = [L(1) 2*LtotLayer+L(2)];
%     
%     switch optOut
%         case 'full'
%             varOut  = varIn;
%         case 'part'
%             if 2*Lbuffer+3*L(2) < Ltot(2) - LAbsLayer
%                 varOut = varIn(:,LAbsLayer+1-L(2):LtotLayer+2*L(2)+Lbuffer,:,:);
%             else
%                 varOut = varIn(:,LAbsLayer/2+1:end-LAbsLayer/2,:,:);
%             end                    
%         case 'guide'
%             varOut  = varIn(:,LAbsLayer+1:LtotLayer+L(2)+Lbuffer,:,:);
%         case 'sample'
%             varOut  = varIn(:,LtotLayer+1:LtotLayer+L(2),:,:);
%         case 'detector'
%             detectorPos = paramLattice.detectorPos;
%             detectorPos = cell2mat(detectorPos);
%             indDetector = sub2ind(Ltot,detectorPos(:,1),detectorPos(:,2));
%             indDetector = sort(indDetector,'ascend');
%             indDetector = unique(indDetector);
%             nDetector = length(indDetector);
%             
%             nNN = paramLattice.nNN;
% 
%             varOut = zeros(1,nDetector,nNN,2);
%             for indD = 1:nDetector
%                 [indi,indj] = ind2sub(Ltot,indDetector(indD));
%                 varOut(1,indD,:,:) = varIn(indi,indj,:,:);
%             end
%         otherwise
%             error('paramSolver.fullOutput: wrong attribute, needs to be among full, sample or detector');
%     end
% 
% end
% 
% function phi = phiWCA(x,xc)
% 
%     xLJ = 2^(1/6);
%     phi = zeros(size(x));
%     tmpx = (x(x<xc)*xLJ/xc);
%     phi(x < xc) = 4*(tmpx.^-12 - tmpx.^-6)+1;
%     
% end

function [Ec, Ep] = energy(distNode,dU,paramLattice)

    %%
    %%%INPUT
    kMat = paramLattice.kMat;
    mMat = paramLattice.mMat;
    m0 = 1/(2*pi)^2;
    a = paramLattice.a;
    LAbsLayer = paramLattice.LAbsLayer;
    LtotLayer = paramLattice.LtotLayer;
    Lbuffer = paramLattice.Lbuffer;
    L = paramLattice.L;

    %%%ENERGY
    Ec = m0*sum(dU.^2,3)/2;
    Ep = 1/2*sum(kMat.*(distNode - a).^2,3)/2;
    Ec(mMat == 0) = 0;
    Ep(mMat == 0) = 0;

    try
        switch paramLattice.removeLead
            case 'left'
                Ec = Ec(:,1:L(2)+Lbuffer);
                Ep = Ep(:,1:L(2)+Lbuffer);
            case 'right'
                Ec = Ec(:,LAbsLayer+1:LtotLayer+L(2));
                Ep = Ep(:,LAbsLayer+1:LtotLayer+L(2));
        end
    catch
        Ec = Ec(:,LAbsLayer+1:LtotLayer+L(2)+Lbuffer);
        Ep = Ep(:,LAbsLayer+1:LtotLayer+L(2)+Lbuffer);
    end


end

















































