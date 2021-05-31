function [scatMat, fSampling] = scatteringMatrixPulseSHengine6(U0,paramLattice,f,energyThreshold,fRange,modeRange,optDisp)


    %%
    %%%%%%%%%%
    %VARIABLES
    %%%%%%%%%%
    L = paramLattice.L;
    L = L(1);
    
    amplitude = 0.001;
    forcing = 'position';
    optBC = 'clamped';
    periodSampling = 50;
    freqSave = 1;
    stopOpt = 'energy';
    nIter = 5000;
    
    paramSolver = genParamSolverRelax(periodSampling,nIter,stopOpt,energyThreshold,freqSave,f,optBC);
    paramSolver.optOut = 'detector';
    paramSolver.optDisp = optDisp;
    paramSolver.maxIter = 500000;
    paramSolver.minIter = 10000;
    
    df = 0.00001;
    fSampling = fRange(1):df:fRange(2);
    nf = length(fSampling);
    
    nmin = ceil(modeRange(1));
    nmax = floor(modeRange(2));
    
    guideBasis = (nmin:nmax)';
    nModesIn = length(guideBasis);
    guideBasis = [guideBasis ones(nModesIn,1); guideBasis 2*ones(nModesIn,1)];
    
    waveVec = waveVecGuideMode(fSampling,paramLattice,modeRange);
    weightMode = sqrt(real(sin(waveVec)));
    nModesOut = size(weightMode,2);
%     weightMode = weightMode(:,1:nModesIn);
    nCycles = 50;
    
    newfig(13); %monitor
    subplot(2,1,1);
    subplot(2,1,2);
    formatfig(13);
    %%%%%%%%%
    %OUTPUT
    %%%%%%%%%
    
    
    %%
    %%%%%%%%%
    %COMPUTE
    %%%%%%%%%
    scatMat = zeros(2*nModesOut,2*nModesIn,nf);
    %%
    disp('COMPUTE');
    
    parfor indBasisIn = 1:2*nModesIn
        
        %% PREPARE
%         fprintf('%d / %d\n',indBasisIn,2*nModes);
        %%%MODE
        modeIdxIn = guideBasis(indBasisIn,:);
        sourceId = modeIdxIn(2);
        modeIdxIn = modeIdxIn(1);
        
        %%%PARAM SOURCE
        sourcePos = paramLattice.detectorPos{sourceId};
        shapeWave = genAmpGuideMode(modeIdxIn,L,'S');

        tau = nCycles/f;
        paramSource = genParamSourceWavePacket(f,amplitude,nCycles,tau,U0,shapeWave,sourcePos,forcing,paramSolver);

        %% SOLVE
        Ui = zeros(size(U0));
        Ui = Ui(:,:,1);
        dUi = Ui;
        dataOut = verletSolverSHRelax(Ui,dUi,U0,paramLattice,paramSolver,paramSource);
        
        if 0 == 1
            %%
            figure(13);
            close(13);
            figure(13);
            U = dataOut.U;
            subplot(2,1,1);
            plot(U(:,L/2));
            subplot(2,1,1);
            plot(U(:,L/2*3));
        end

        %% PROJECTION
        
        %%PROJECTION
        [shapeSignalIn, shapeSignalOut] = shapeSignalsFFTSHwavePacket2(dataOut,tau,fSampling,paramLattice,paramSolver);
        [coeffIn,  shapeRecoverIn]  = modeProjectionSH2(shapeSignalIn,nmax);
        [coeffOut, shapeRecoverOut] = modeProjectionSH2(shapeSignalOut,nmax);

        
        if 0 == 1
            %%
            figure(13);
            close(13);
            figure(13);
            subplot(2,1,1);
            hold all;
            plot(abs(shapeSignalIn(:,1,(nf+1)/2)));
            plot(abs(shapeSignalOut(:,1,(nf+1)/2)));
            plot(abs(shapeSignalOut(:,2,(nf+1)/2)));
            legend('in','out left','out right');
            subplot(2,1,2);
            hold all;
            plot(abs(shapeRecoverIn(:,1,(nf+1)/2)));
            plot(abs(shapeRecoverOut(:,1,(nf+1)/2)));
            plot(abs(shapeRecoverOut(:,2,(nf+1)/2)));
            
            newfig(14);
            tmpIn  = coeffIn(:,:,(nf+1)/2);
            tmpOut = coeffOut(:,:,(nf+1)/2);
            plot(abs(tmpIn(:)));
            plot(abs(tmpOut(:)));
            legend('in','out')
        end
        
        %%
        coeffIn  = coeffIn(nmin:nmax,:,:);
        coeffOut = coeffOut(nmin:nmax,:,:);
        
        
        %% POST PROCESSING
        
        cIn = permute(coeffIn,[3 1 2]);
        cOut = permute(coeffOut,[3 1 2]);
        
        %%WEIGHT
        cIn  = weightMode.*cIn;
        cOut = weightMode.*cOut;
        
        
        %% COLLECT
        tmpScatMat = zeros(2*nModesOut,nf);
        for indf = 1:nf
            tmpcOut = cOut(indf,:,:);
            tmpcIn  = cIn(indf,:,:);
            tmpInd  = mod(indBasisIn-1,nModesIn)+1;
            tmp = tmpcOut(:)/tmpcIn(1,tmpInd,sourceId);
            tmpScatMat(:,indf) = tmp;
        end
        
        %DISPLAY
        [~,indf0] = min(abs(fSampling - f));
        tmp = tmpScatMat(:,indf0);
        
        totIter = length(dataOut.E);
        stopVar = dataOut.stopVar;
        fprintf('%d / %d - sum coeff^2 = %1.5f\n%d / %d - nIter = %d\tRemaining Energy = %1.2e\n',indBasisIn,2*nModesIn,sum(abs(tmp(:)).^2),indBasisIn,2*nModesIn,totIter,stopVar);
        
        %%%%%
        scatMat(:,indBasisIn,:)  = reshape(tmpScatMat,size(scatMat(:,indBasisIn,:)));
        %%%%%

        if 0 == 1
            %%
            figure(13);
            sumScatMat = sum(abs(tmpScatMat).^2,1);
            subplot(2,1,1);
            cla;
            plot(fSampling, sumScatMat);
            ylim([0 2]);
            subplot(2,1,2);
            
            cla;
            plot(fSampling, sumScatMat);
            ylim([0.95 1.05]);
            formatfig(13);
            drawnow;
        end
        
    end
    


        
end

function [shapeSignalIn, shapeSignalOut] = shapeSignalsFFTSHwavePacket2(data,timeThreshold,fSampling,paramLattice,paramSolver)
    

    %%PARAMETERS
    UvsTime = data.U;
    time = data.time;

    L = paramLattice.L;
    LtotLayer = paramLattice.LtotLayer;
    Lbuffer = paramLattice.Lbuffer;
    LAbsLayer = paramLattice.LAbsLayer;
    Ltot = [L(1) L(2)+2*LtotLayer];
    detectorPos = paramLattice.detectorPos;
    

    
    df = fSampling(2) - fSampling(1);
    nf = length(fSampling);
    dt = time(2) - time(1);
    nt = round(1/(df*dt));
    if nt < length(time)
        error('Input frequency sampling is not dense enought to cover FFT padding');
    end
    
    %% FORMAT DATA
    
    boolOut = time > timeThreshold;
    boolIn = time <= timeThreshold;
    timeOut = time(boolOut);
    timeIn = time(boolIn);
    
    switch paramSolver.optOut
        case 'detector'
            indPos = detectorPos2ind(detectorPos,Ltot);
            indL = indPos{1};
            indR = indPos{2};

            UoutL = UvsTime(boolOut,indL,:);
            UoutR = UvsTime(boolOut,indR,:);
            UinL = UvsTime(boolIn,indL,:);
            UinR = UvsTime(boolIn,indR,:);
        case 'sample'
            error('detectors out of recorded ROI. OptOut needs to be either ''detector'', ''guide'', ''part'' or ''full''');
        case 'guide'
            UoutL = UvsTime(boolOut,:,1,:);
            UoutR = UvsTime(boolOut,:,end,:);
            UinL = UvsTime(boolIn,:,1,:);
            UinR = UvsTime(boolIn,:,end,:);
        case 'part'
            if 2*Lbuffer+3*L(2) < Ltot(2) - LAbsLayer
                UoutL = UvsTime(boolOut,:,1+L(2),:);
                UoutR = UvsTime(boolOut,:,end-L(2),:);
                UinL = UvsTime(boolIn,:,1+L(2),:);
                UinR = UvsTime(boolIn,:,end-L(2),:);
            else
                UoutL = UvsTime(boolOut,:,LAbsLayer/2+1,:);
                UoutR = UvsTime(boolOut,:,end-LAbsLayer/2,:);
                UinL = UvsTime(boolIn,:,LAbsLayer/2+1,:);
                UinR = UvsTime(boolIn,:,end-LAbsLayer/2,:);
            end    
        case 'full'
                UoutL = UvsTime(boolOut,:,LAbsLayer/2+1,:);
                UoutR = UvsTime(boolOut,:,end-LAbsLayer/2,:);
                UinL = UvsTime(boolIn,:,LAbsLayer/2+1,:);
                UinR = UvsTime(boolIn,:,end-LAbsLayer/2,:);
        otherwise
            error('CACA paramSolver.optOut has wrong attribute');
    end
    UoutL = squeeze(UoutL);
    UoutR = squeeze(UoutR);
    UinL = squeeze(UinL);
    UinR = squeeze(UinR);
        
    %% FOURIER TRANSFORM
    optNorm = 'energy';
    fftUoutL = zeros(nf,L(1));
    fftUoutR = zeros(nf,L(1));
    fftUinL  = zeros(nf,L(1));
    fftUinR  = zeros(nf,L(1));
    for indL = 1:L(1)

        [tmpfFFT,tmpFFT] = fourierTransform(timeOut,UoutL(:,indL),nt,optNorm);
        fftUoutL(:,indL) = interp1(tmpfFFT,tmpFFT,fSampling,'spline');
        
        [tmpfFFT,tmpFFT] = fourierTransform(timeOut,UoutR(:,indL),nt,optNorm);
        fftUoutR(:,indL) = interp1(tmpfFFT,tmpFFT,fSampling,'spline');
        
        [tmpfFFT,tmpFFT] = fourierTransform(timeOut,UinL(:,indL),nt,optNorm);
        fftUinL(:,indL) = interp1(tmpfFFT,tmpFFT,fSampling,'spline');
        
        [tmpfFFT,tmpFFT] = fourierTransform(timeOut,UinR(:,indL),nt,optNorm);
        fftUinR(:,indL) = interp1(tmpfFFT,tmpFFT,fSampling,'spline');
        
    end
    %%

    
    shapeSignalIn = zeros(L(1),2,nf);
    shapeSignalOut = zeros(L(1),2,nf);
    shapeSignalIn(:,1,:) = reshape(fftUinL',[L(1) 1 nf]);
    shapeSignalIn(:,2,:) = reshape(fftUinR',[L(1) 1 nf]);
    shapeSignalOut(:,1,:) = reshape(fftUoutL',[L(1) 1 nf]);
    shapeSignalOut(:,2,:) = reshape(fftUoutR',[L(1) 1 nf]);

    
    
end


function [projCoeff, shapeRecover] = modeProjectionSH2(shapeIn,nModes)

    L = size(shapeIn,1);
    nf = size(shapeIn,3);
    
    projCoeff = zeros(nModes,2,nf);
    shapeRecover = 0*shapeIn;
    
    for indf = 1:nf
        for indLR = 1:2
            for indMode = 1:nModes
                waveShape = genAmpGuideMode(indMode,L(1),'S');
                projCoeff(indMode,indLR,indf) = mean(sum(shapeIn(:,indLR,indf).*waveShape,1))./sum(waveShape.^2);
            end
        end
    end
    
    
    %%
    for indf = 1:nf
        shapeOutL = zeros(L,1);
        shapeOutR = zeros(L,1);
        for indMode = 1:nModes
            waveShape = genAmpGuideMode(indMode,L(1),'S');
            shapeOutL = shapeOutL + projCoeff(indMode,1,indf)*waveShape;
            shapeOutR = shapeOutR + projCoeff(indMode,2,indf)*waveShape;
        end
        
        shapeRecover(:,1,indf) = shapeOutL;
        shapeRecover(:,2,indf) = shapeOutR;
    end


   
    
    
    
    
end






























