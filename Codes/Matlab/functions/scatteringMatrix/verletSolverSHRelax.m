function dataOut = verletSolverSHRelax(Ui,dUi,U0,paramLattice,paramSolver,paramSource)
    % Propagates a wave in a 2D masses and springs percolated network

    stopOpt = paramSolver.stopOpt;
    stopThreshold = paramSolver.stopThreshold;
    nIter = paramSolver.nIter;
    optDisp = paramSolver.optDisp;
    paramSolver.optDisp = 0;
    minIter = paramSolver.minIter;
    maxIter = paramSolver.maxIter;
    optEflux = paramSolver.optEflux;
    
    %%Monitor Energy
    [Ec, Ep] = energyMapSH(Ui,dUi,U0,paramLattice);
    E0 = [sum(Ec(:)) sum(Ep(:))];
    
    indt = 0;
    varMonitor = 1;
    switch stopOpt
        case 'energy'
            maxE = 0;
%         case 'flux'
%             maxFlux = 0;
    end
    %Source Signal
    USources = paramSource.USources;
    
    %%VAR OUT
    UvsTime = [];
    dUvsTime = [];
    EvsTime = [];
    forceSourceVStime = [];
    Eflux = [];

    %% LOOP
    while (varMonitor > stopThreshold || indt < minIter) && ~(indt >= maxIter)
        
        if 0 == 1
            %%
            figure(1);
            clf;
            subplot(2,1,1);
            imagesc(Ui);
            subplot(2,1,2);
            imagesc(dUi);
            pause;
            
        end
        
        %% SOLVE
        if indt <= size(USources,1) - nIter
            paramSource.USources = USources(indt+1:indt+nIter,:);
        elseif indt <= size(USources,1)
            paramSource.USources = USources(indt+1:end,:);
        else
            paramSource.USources = [];
        end
        
        
        tmpDataOut = verletSolverSH(Ui,dUi,U0,paramLattice,paramSolver,paramSource);

        %%
        Uf = tmpDataOut.Ulast;
        dUf = tmpDataOut.dUlast;
        E = tmpDataOut.E;
        Ef = E(end,:);

        Uf = squeeze(Uf);
        dUf = squeeze(dUf);
        Ef = Ef - E0;
        Ef = sum(Ef);
        
        %% MONITOR

        indt = indt + nIter;
        switch stopOpt
            case 'disp'
                varMonitor = sqrt(max(max(sum((Uf-Ui).^2,3))));
                if optDisp
                    fprintf('%d - Max Disp %2.4f%% / %2.3f%%\n',indt,100*varMonitor,100*stopThreshold);
                end
            case 'energy'
                
                if indt > size(USources,1)
                    if Ef > maxE 
                        maxE = Ef;
                    end
                    varMonitor = Ef/maxE;
                end
                
                if optDisp
                    fprintf('%d - Energy %1.2e - Remain %2.2e / %2.2e\n',indt,Ef,varMonitor,stopThreshold);
                end
        end
%         toc;
        
        if 0 == 1
            %%

            indDir = 1;
            optData = 'position';
            paramVideo = genParamVideo(paramSolver);
            paramVideo.freqFrame = 10;
            paramVideo.UScale = 0.2;
            dispVideo(2,tmpDataOut,indDir,optData,paramLattice,paramVideo,paramSource);
            
        end
        
        
        %% VAR OUT
        
        tmpU  = tmpDataOut.U;
        tmpdU = tmpDataOut.dU;
        tmpE  = tmpDataOut.E;
%         tmpF  = tmpDataOut.forceSource;
        
        UvsTime  = cat(1,UvsTime,tmpU);
        dUvsTime = cat(1,dUvsTime,tmpdU);
        EvsTime  = cat(1,EvsTime,tmpE);
%         forceSourceVStime = cat(1,forceSourceVStime,tmpF);
        
        if optEflux
            tmpflux = tmpDataOut.flux;
            Eflux = cat(1,Eflux,tmpflux);
        end
        %% NEXT
        
        Ui = Uf;
        dUi = dUf;
        
        
    end

    freqSave = paramSolver.freqSave;
    dataOut.U = UvsTime;
    dataOut.dU = dUvsTime;
    dataOut.E = EvsTime;
    dataOut.time = (1:freqSave:indt)*paramSolver.dt;
    dataOut.stopVar = varMonitor;
    dataOut.maxE = maxE;
%     dataOut.forceSource = forceSourceVStime;
    if optEflux
        dataOut.flux = Eflux;
    end
    
    if indt >= maxIter
        fprintf('Simulation run reached maximum iteration limit\n');
    end
    
end




















