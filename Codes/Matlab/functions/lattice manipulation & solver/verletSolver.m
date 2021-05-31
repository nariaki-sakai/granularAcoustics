function dataOut = verletSolver(U0,dU0,paramLattice,paramSolver,paramSD)


    
    %%%%%%%%%%%
    %PARAMETERS
    %%%%%%%%%%%
    m0 = 1/(2*pi)^2;
    LAbsLayer = paramLattice.LAbsLayer;
    LtotLayer = paramLattice.LtotLayer;
    Lbuffer = LtotLayer - LAbsLayer;
    L = paramLattice.L;
    
    nNN = paramLattice.nNN;
    
    % paramSolver
    dt = paramSolver.dt;
    nIter = paramSolver.nIter;
    freqSave = paramSolver.freqSave;
    nt = floor(nIter/freqSave);
    optDisp = paramSolver.optDisp;    
    
    
    %%%%%%%
    %OUTPUT
    %%%%%%%
    UvsTime = formatVarOut3(U0,paramLattice,paramSolver.optOut);
    UvsTime = zeros([nt size(UvsTime)]);
    sizeU = size(UvsTime);
        
    dUvsTime = zeros(sizeU);
    if paramSolver.optEflux
        Eflux = nan([sizeU(1:3) nNN 2]);
    end
    
    EvsTime = zeros(nt,2);
    time = (1:nt)*freqSave*dt;
    
    %%%%%%%%%%%
    %INITIALIZE
    %%%%%%%%%%%
    Ux = U0(:,:,1);
    Uy = U0(:,:,2);
    dUx = dU0(:,:,1);
    dUy = dU0(:,:,2);

    Unx  = neighbourVal(Ux);
    Uny  = neighbourVal(Uy);
    dUnx = neighbourVal(dUx);
    dUny = neighbourVal(dUy);
    [FxC, FyC] = forceConserv(Ux,Uy,Unx,Uny,paramLattice);
    [FxD, FyD] = forceDissip(dUx,dUy,dUnx,dUny,paramLattice);
    d2Ux1 = (FxC + FxD)/m0;
    d2Uy1 = (FyC + FyD)/m0;
    
    U = zeros([size(Ux) 2]);
    dU = U;
    
    %%
    for indt = 1:nIter

        if optDisp
            dispIndForLoop(indt,nIter,10^(round(log10(nIter/10))));
        end
        %% VERLET
        %POSITION
        Ux = Ux + dUx*dt + d2Ux1*dt^2/2;
        Uy = Uy + dUy*dt + d2Uy1*dt^2/2;

        %VELOCITY PASS 0
        tmpdUx = dUx + d2Ux1*dt;
        tmpdUy = dUy + d2Uy1*dt;
        
        if strcmp(paramSD.forcing,'position')
            [Ux,Uy,tmpdUx,tmpdUy] = fixSources(Ux,Uy,tmpdUx,tmpdUy,paramSD,paramSolver,indt);
        end
        [Ux,Uy,tmpdUx,tmpdUy] = fixBoundaryConditions(Ux,Uy,tmpdUx,tmpdUy,paramSolver,paramLattice);
        
        %ACCELERATION AT T+1 PASS 1
        Unx  = neighbourVal(Ux);
        Uny  = neighbourVal(Uy);
        dUnx = neighbourVal(tmpdUx);
        dUny = neighbourVal(tmpdUy);
        
        [FxC, FyC, distNode, FxCAll, FyCAll] = forceConserv(Ux,Uy,Unx,Uny,paramLattice);
        [FxD, FyD] = forceDissip(tmpdUx,tmpdUy,dUnx,dUny,paramLattice);
        d2Ux2 = (FxC + FxD)/m0;
        d2Uy2 = (FyC + FyD)/m0;
        
        if strcmp(paramSD.forcing,'force')
            [d2Ux2,d2Uy2] = addForcing(d2Ux2,d2Uy2,paramSD,indt);
        end

        %VELOCITY PASS 1
        tmpdUx = dUx + (d2Ux1+d2Ux2)*dt/2;
        tmpdUy = dUy + (d2Uy1+d2Uy2)*dt/2;
        
        if strcmp(paramSD.forcing,'position')
            [Ux,Uy,tmpdUx,tmpdUy] = fixSources(Ux,Uy,tmpdUx,tmpdUy,paramSD,paramSolver,indt);
        end
        [Ux,Uy,tmpdUx,tmpdUy] = fixBoundaryConditions(Ux,Uy,tmpdUx,tmpdUy,paramSolver,paramLattice);
        
        %ACCELERATION AT T+1 PASS 2
        dUnx = neighbourVal(tmpdUx);
        dUny = neighbourVal(tmpdUy);
        [FxD, FyD, FxFAll, FyFAll, FxV, FyV] = forceDissip(tmpdUx,tmpdUy,dUnx,dUny,paramLattice);
        d2Ux2 = (FxC + FxD)/m0;
        d2Uy2 = (FyC + FyD)/m0;
        
        if strcmp(paramSD.forcing,'force')
            [d2Ux2,d2Uy2] = addForcing(d2Ux2,d2Uy2,paramSD,indt);
        end
        
        %VELOCITY PASS 2
        dUx = dUx + (d2Ux1+d2Ux2)*dt/2;
        dUy = dUy + (d2Uy1+d2Uy2)*dt/2;
        
        if strcmp(paramSD.forcing,'position')
            [Ux,Uy,dUx,dUy] = fixSources(Ux,Uy,dUx,dUy,paramSD,paramSolver,indt);
        end
        [Ux,Uy,dUx,dUy] = fixBoundaryConditions(Ux,Uy,dUx,dUy,paramSolver,paramLattice);
        
        %%NEXT
        d2Ux1 = d2Ux2;
        d2Uy1 = d2Uy2;
        
        %% Save
        if mod(indt,freqSave) == 0

            U(:,:,1) = Ux;
            U(:,:,2) = Uy;
            dU(:,:,1) = dUx;
            dU(:,:,2) = dUy;
            dUn = zeros([size(dUnx) 2]);
            dUn(:,:,:,1) = neighbourVal(dUx);
            dUn(:,:,:,2) = neighbourVal(dUy);
            
            if paramSolver.optEflux
                flux = energyFlux(dU,dUn,FxCAll,FyCAll,FxFAll,FyFAll,FxV,FyV,paramLattice);
                Eflux(floor(indt/freqSave),:,:,:,:)  = formatVarOut4(flux,paramLattice,paramSolver.optOut);
            end
            [Ec, Ep] = energyMap(distNode,dU,paramLattice);
            Ec = Ec(:,LAbsLayer+1:LtotLayer+L(2)+Lbuffer);
            Ep = Ep(:,LAbsLayer+1:LtotLayer+L(2)+Lbuffer);
            
            UvsTime(floor(indt/freqSave),:,:,:)  = formatVarOut3(U,paramLattice,paramSolver.optOut);
            dUvsTime(floor(indt/freqSave),:,:,:) = formatVarOut3(dU,paramLattice,paramSolver.optOut);
            EvsTime(floor(indt/freqSave),1) = sum(Ec(:));
            EvsTime(floor(indt/freqSave),2) = sum(Ep(:));

        end
        
    end
    
    UvsTime = squeeze(UvsTime);
    dUvsTime = squeeze(dUvsTime);
    
    dataOut.U = UvsTime;
    dataOut.dU = dUvsTime;
    dataOut.time = time;
    
    U0 = formatVarOut3(U0,paramLattice,paramSolver.optOut);
    U0 = squeeze(U0);
    U0 = reshape(U0,[1 size(U0)]);
    dataOut.U0 = U0;
    
    dataOut.E = EvsTime;
    if paramSolver.optEflux
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

end

function [distNodeX, distNodeY, distNode] = distanceNode(Ux,Uy,Unx,Uny)

    distNodeX = zeros(size(Unx));
    distNodeY = zeros(size(Uny));
    distNode  = zeros(size(Unx));
    
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
        distNode(:,:,indNN) = sqrt(distNodeX(:,:,indNN).^2 + distNodeY(:,:,indNN).^2);
    end
    for indNN = 5:8
        distNodeX(:,:,indNN) = -circshift(distNodeX(:,:,9-indNN),vecLattice(9-indNN,:));
        distNodeY(:,:,indNN) = -circshift(distNodeY(:,:,9-indNN),vecLattice(9-indNN,:));
        distNode(:,:,indNN)  =  circshift(distNode(:,:,9-indNN),vecLattice(9-indNN,:));
    end
% 
%         %%Springs lengths
%     [ lDLx,lDLy ,lDx,lDy,lDRx,lDRy,lRx,lRy,lURx,lURy,lUx,lUy,lULx,lULy,lLx,lLy,lDL,lD,lDR,lR,lUR,lU,lUL,lL] = lengths_vectors(Ux,Uy,1);
%     distNodeX(:,:,1) = lDLx;
%     distNodeX(:,:,2) = lLx;
%     distNodeX(:,:,3) = lULx;
%     distNodeX(:,:,4) = lDx;
%     distNodeX(:,:,5) = lUx;
%     distNodeX(:,:,6) = lDRx;
%     distNodeX(:,:,7) = lRx;
%     distNodeX(:,:,8) = lURx;
%     
%     distNodeY(:,:,1) = lDLy;
%     distNodeY(:,:,2) = lLy;
%     distNodeY(:,:,3) = lULy;
%     distNodeY(:,:,4) = lDy;
%     distNodeY(:,:,5) = lUy;
%     distNodeY(:,:,6) = lDRy;
%     distNodeY(:,:,7) = lRy;
%     distNodeY(:,:,8) = lURy;
%     
%     distNode(:,:,1) = lDL;
%     distNode(:,:,2) = lL;
%     distNode(:,:,3) = lUL;
%     distNode(:,:,4) = lD;
%     distNode(:,:,5) = lU;
%     distNode(:,:,6) = lDR;
%     distNode(:,:,7) = lR;
%     distNode(:,:,8) = lUR;
        
end

function [FxC, FyC, distNode, FxCAll, FyCAll] = forceConserv(Ux,Uy,Unx,Uny,paramLattice)

    %%PARAMETERS
%     strain = paramLattice.strain;
    kMat = paramLattice.kMat;
    mMat = paramLattice.mMat;
    a = paramLattice.a;

    %%DIST NODES
    [distNodeX, distNodeY, distNode] = distanceNode(Ux,Uy,Unx,Uny);
        
    %%X
    FxCAll = kMat.*(1-a./distNode).*distNodeX;
    FxCAll(distNode==0) = 0;
    FxC = sum(FxCAll,3);
    FxC(mMat == 0) = 0;
    dispNaNInf(FxC);
    
    %%Y
    FyCAll = kMat.*(1-a./distNode).*distNodeY;
    FyCAll(distNode==0) = 0;
    FyC = sum(FyCAll,3);
    FyC(mMat == 0) = 0;
    dispNaNInf(FyC);
        
end

function [FxD, FyD, FxFAll, FyFAll, FxV, FyV] = forceDissip(dUx,dUy,dUnx,dUny,paramLattice)

    %%PARAMETERS
    viscousMat = paramLattice.viscous;
    frictionMat = paramLattice.friction;
    m0 = 1/(2*pi)^2;
    mMat = paramLattice.mMat;

    %%VISCOUS
    FxV = - viscousMat .* dUx;
    FyV = - viscousMat .* dUy;
    FxV(mMat == 0) = 0;
    FyV(mMat == 0) = 0;
    
    %%% FRICTION X
    dUx0 = reshape(dUx,[size(dUx) 1]);
    relVelX = dUnx - dUx0;
    FxFAll = frictionMat.*relVelX;
    FxF = sum(FxFAll,3);
    FxF(mMat == 0) = 0;
    
    %%% FRICTION Y
    dUy0 = reshape(dUy,[size(dUy) 1]);
    relVelY = dUny - dUy0;
    FyFAll = m0*frictionMat.*relVelY;
    FyF = sum(FyFAll,3);
    FyF(mMat == 0) = 0;

    %%% CHECK AND OUTPUT
    dispNaNInf(FxF);
    dispNaNInf(FxV);
    dispNaNInf(FyF);
    dispNaNInf(FyV);
    FxD = FxF + FxV;
    FyD = FyF + FyV;

end

function [Ux,Uy,dUx,dUy] = fixSources(Ux,Uy,dUx,dUy,paramSD,paramSolver,indt)
    
    dt = paramSolver.dt;
    try 
        USources = paramSD.USources;
        sourcePos = paramSD.sourcePos;
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

            Ux(indi,indj)  =  USources(indt,indSource,1);
            Uy(indi,indj)  =  USources(indt,indSource,2);
            if indt > 1 && indt < ntSource
                dUx(indi,indj)  =  (USources(indt+1,indSource,1)-USources(indt-1,indSource,1))/(2*dt);
                dUy(indi,indj)  =  (USources(indt+1,indSource,2)-USources(indt-1,indSource,2))/(2*dt);
            elseif indt == 1
                dUx(indi,indj)  =  (USources(indt+1,indSource,1)-USources(indt,indSource,1))/dt;
                dUy(indi,indj)  =  (USources(indt+1,indSource,2)-USources(indt,indSource,2))/dt;
            else
                dUx(indi,indj)  =  (USources(indt,indSource,1)-USources(indt-1,indSource,1))/dt;
                dUy(indi,indj)  =  (USources(indt,indSource,2)-USources(indt-1,indSource,2))/dt;
            end
        end
    end
    
end

function [Ux,Uy,dUx,dUy] = fixBoundaryConditions(Ux,Uy,dUx,dUy,paramSolver,paramLattice)
    %%

    L = paramLattice.L;
    strain = paramLattice.strain;
    mMat = paramLattice.mMat;
    
    LtotLayer = paramLattice.LtotLayer;
    Ltot = [L(1) L(2)+2*LtotLayer];
    x = 1:Ltot(2);
    y = 1:Ltot(1);
    wallDispX = x*strain;
    wallDispY = y*strain;

    fixedBC = paramSolver.fixedBC;
    if fixedBC(1) == 1
        boolMass = mMat(:,1) > 0;
        Ux(boolMass,1)  = wallDispX(1);
        Uy(boolMass,1)  = wallDispY(boolMass)'; 
        dUx(:,1) = 0;
        dUy(:,1) = 0; 
    end
    if fixedBC(2) == 1
        boolMass = mMat(1,:) > 0;
        Ux(1,boolMass)  = wallDispX(boolMass);
        Uy(1,boolMass)  = wallDispY(1);
        dUx(1,boolMass) = 0;
        dUy(1,boolMass) = 0;    
    end
    if fixedBC(3) == 1
        boolMass = mMat(end,:) > 0;
        Ux(end,boolMass)  = wallDispX(boolMass);
        Uy(end,boolMass)  = wallDispY(end); 
        dUx(end,boolMass) = 0;
        dUy(end,boolMass) = 0;
    end
    if fixedBC(4) == 1
        boolMass = mMat(:,end) > 0;
        Ux(boolMass,end)  = wallDispX(end);
        Uy(boolMass,end)  = wallDispY(boolMass)';
        dUx(:,end) = 0;
        dUy(:,end) = 0;
    end

end
% 
% function [Ux,Uy,dUx,dUy] = fixBoundaryConditions(Ux,Uy,dUx,dUy,paramSolver,paramLattice)
%     
%     wallPos = paramSolver.fixedBC;
%     
%     Ux(:,1)  = wallPos(1,1);
%     Uy(1,:)  = wallPos(2,2);
%     Uy(end,:)  = wallPos(3,2);
%     Ux(:,end)  = wallPos(4,1);        
%     dUx(:,1) = 0;
%     dUy(1,:) = 0;
%     dUy(end,:) = 0;  
%     dUx(:,end) = 0;
%     
% end

function [d2Ux,d2Uy] = addForcing(d2Ux,d2Uy,paramSD,indt)

    m0 = 1/(2*pi)^2;
    try 
        USources = paramSD.USources;
        sourcePos = paramSD.sourcePos;
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

            d2Ux(indi,indj) = d2Ux(indi,indj) + USources(indt,indSource,1)/m0;
            d2Uy(indi,indj) = d2Uy(indi,indj) + USources(indt,indSource,2)/m0;
        end
    end

end

function Eflux = energyFlux(dU,dUn,FxCAll,FyCAll,FxFAll,FyFAll,FxV,FyV,paramLattice)

    %%
    %%%INPUT
    Ltot = size(dU);
    try
        nNN = paramLattice.nNN;
    catch
        nNN = 4;
    end

    %%%OUTPUT
    Eflux = zeros(Ltot(1),Ltot(2),nNN,2);

    %%%FLUX
    if nNN == 4
        
        for ind = 1:4
            indNN = 2*ind;
            if indNN > 5
                indNN = indNN - 1;
            end
            Eflux(:,:,ind,1) = FxCAll(:,:,indNN).*(dUn(:,:,indNN,1)+dU(:,:,1))/2;% + ...
                               %FxFAll(:,:,indNN).*dU(:,:,1) + ...
                               %FxV(:,:,1).*dU(:,:,1);
            Eflux(:,:,ind,2) = FyCAll(:,:,indNN).*(dUn(:,:,indNN,2)+dU(:,:,2))/2;% + ...
                               %FyFAll(:,:,indNN).*dU(:,:,2) + ...
                               %FyV(:,:).*dU(:,:,2);
        end
        
    elseif nNN == 8
        for indNN = 1:8
            Eflux(:,:,indNN,1) = FxCAll(:,:,indNN).*(dUn(:,:,indNN,1)+dU(:,:,1))/2;% + ...
%                                  FxFAll(:,:,indNN).*dU(:,:,1) + ...
%                                  FxV(:,:,1).*dU(:,:,1);
            Eflux(:,:,indNN,2) = FyCAll(:,:,indNN).*(dUn(:,:,indNN,2)+dU(:,:,2))/2;% + ...
%                                  FyFAll(:,:,indNN).*dU(:,:,2) + ...
%                                  FyV(:,:).*dU(:,:,2);
        end
        
    else
        error('number of neighbours should be 4 or 8');
    end

end

function varOut = formatVarOut3(varIn,paramLattice,optOut)

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
            
            varOut = zeros(1,nDetector,2);
            for indD = 1:nDetector
                [indi,indj] = ind2sub(Ltot,indDetector(indD));
                varOut(1,indD,:) = varIn(indi,indj,:);
            end
        otherwise
            error('paramSolver.fullOutput: wrong attribute, needs to be among full, sample or detector');
    end

end

function varOut = formatVarOut4(varIn,paramLattice,optOut)

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
                varOut = varIn(:,LAbsLayer+1-L(2):LtotLayer+2*L(2)+Lbuffer,:,:);
            else
                varOut = varIn(:,LAbsLayer/2+1:end-LAbsLayer/2,:,:);
            end                    
        case 'guide'
            varOut  = varIn(:,LAbsLayer+1:LtotLayer+L(2)+Lbuffer,:,:);
        case 'sample'
            varOut  = varIn(:,LtotLayer+1:LtotLayer+L(2),:,:);
        case 'detector'
            detectorPos = paramLattice.detectorPos;
            detectorPos = cell2mat(detectorPos);
            indDetector = sub2ind(Ltot,detectorPos(:,1),detectorPos(:,2));
            indDetector = sort(indDetector,'ascend');
            indDetector = unique(indDetector);
            nDetector = length(indDetector);
            
            nNN = paramLattice.nNN;

            varOut = zeros(1,nDetector,nNN,2);
            for indD = 1:nDetector
                [indi,indj] = ind2sub(Ltot,indDetector(indD));
                varOut(1,indD,:,:) = varIn(indi,indj,:,:);
            end
        otherwise
            error('paramSolver.fullOutput: wrong attribute, needs to be among full, sample or detector');
    end

end

function [Ec, Ep] = energyMap(distNode,dU,paramLattice)

    %%
    %%%INPUT
    kMat = paramLattice.kMat;
    mMat = paramLattice.mMat;
    m0 = 1/(2*pi)^2;
    a = paramLattice.a;

    %%%ENERGY
    Ec = m0*sum(dU.^2,3)/2;
    Ep = 1/2*sum(kMat.*(distNode - a).^2,3)/2;
    Ec(mMat == 0) = 0;
    Ep(mMat == 0) = 0;


end



















































