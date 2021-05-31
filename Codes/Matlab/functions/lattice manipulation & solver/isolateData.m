function [dataOut, paramLatticeOut] = isolateData(dataIn,paramLatticeIn,paramSolver,optOut)

    optIn = paramSolver.optOut;
    L = paramLatticeIn.L;
    Lbuffer = paramLatticeIn.Lbuffer;
    LtotLayer = paramLatticeIn.LtotLayer;
    LAbsLayer = paramLatticeIn.LAbsLayer;
    mMat = paramLatticeIn.mMat;
    kMat = paramLatticeIn.kMat;
    
    Uin = dataIn.U;
    dUin = dataIn.dU;
    if paramSolver.optEflux
        fluxIn = dataIn.flux;
    end
    Ulast = dataIn.Ulast;
    Ltot = size(Ulast);
    nIter = size(Uin,1);
    
    paramLatticeOut = paramLatticeIn;
    
    switch optIn
        case 'full'
        case 'part'
            if 2*Lbuffer+3*L(2) < Ltot(2) - LAbsLayer
                switch optOut
                    case 'part'
                        U  = Uin;
                        dU  = dUin;
                        if paramSolver.optEflux
                            flux  = fluxIn;
                        end
                    case 'guide'
                        U  = Uin(:,:,L(2)+1:end-L(2));
                        dU  = dUin(:,:,L(2)+1:end-L(2));
                        if paramSolver.optEflux
                            flux  = fluxIn(:,:,L(2)+1:end-L(2),:);
                        end
                    case 'sample'
                        U  = Uin(:,:,L(2)+Lbuffer+1:end-L(2)-Lbuffer);
                        dU  = dUin(:,:,L(2)+1+Lbuffer:end-L(2)-Lbuffer);
                        if paramSolver.optEflux
                            flux  = fluxIn(:,:,L(2)+1+Lbuffer:end-L(2)-Lbuffer,:);
                        end
                    otherwise
                        error('optOut has wrong attribute');
                end
            else
                switch optOut
                    case 'part'
                        U  = Uin;
                        dU  = dUin;
                        if paramSolver.optEflux
                            flux  = fluxIn;
                        end
                    case 'guide'
                        U  = Uin(:,:,LAbsLayer/2+1:end-LAbsLayer/2);
                        dU = dUin(:,:,LAbsLayer/2+1:end-LAbsLayer/2);
                        if paramSolver.optEflux
                            flux = fluxIn(:,:,LAbsLayer/2+1:end-LAbsLayer/2,:);
                        end
                    case 'sample'
                        U  = Uin(:,:,LAbsLayer/2+1:end-LAbsLayer/2);
                        dU = dUin(:,:,LAbsLayer/2+1:end-LAbsLayer/2);
                        if paramSolver.optEflux
                            flux = fluxIn(:,:,LAbsLayer/2+1:end-LAbsLayer/2,:);
                        end
                otherwise
                        error('optOut has wrong attribute');
                end

            end 
        case 'guide'
            switch optOut
                case 'sample'
                    U  = Uin(:,:,Lbuffer+1:end-Lbuffer);
                    dU  = dUin(:,:,1+Lbuffer:end-Lbuffer);
                    if paramSolver.optEflux
                        flux  = fluxIn(:,:,1+Lbuffer:end-Lbuffer,:);
                    end
                case 'guide'
                    U  = Uin;
                    dU  = dUin;
                    if paramSolver.optEflux
                        flux  = fluxIn;
                    end
                otherwise
                    error('optOut has wrong attribute');
            end
        otherwise
            error('paramSolver.optOut has wrong attribute');
    end
    
    
    dataOut = dataIn;
    dataOut.U = U;
    dataOut.dU = dU;
    if paramSolver.optEflux
        dataOut.flux = flux;
    end
    
    switch optOut
        case 'part'
            if 2*Lbuffer+3*L(2) < Ltot(2) - LAbsLayer
                mMat  = mMat(:,LAbsLayer-L(2)+1:end-LAbsLayer+L(2));
                kMat  = kMat(:,LAbsLayer-L(2)+1:end-LAbsLayer+L(2),:);
            else
                mMat  = mMat(:,LAbsLayer/2+1:end-LAbsLayer/2);
                kMat  = kMat(:,LAbsLayer/2+1:end-LAbsLayer/2,:);
            end 
        case 'guide'
            mMat = mMat(:,LAbsLayer+1:end-LAbsLayer);
            kMat = kMat(:,LAbsLayer+1:end-LAbsLayer,:);
        case 'sample'
            mMat = mMat(:,LtotLayer+1:end-LtotLayer);
            kMat = kMat(:,LtotLayer+1:end-LtotLayer,:);
        otherwise
            error('optOut has wrong attribute');
    end
    paramLatticeOut.mMat = mMat;
    paramLatticeOut.kMat = kMat;
    

end

