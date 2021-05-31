function paramSource = genParamSourceMonoChrome(f,sourceAmp,nCycles,tau,U0,coeffWave,sourcePos,forcing,paramSolver)

    %%
    L = size(U0);
    dt = paramSolver.dt;
    
    nSources = size(sourcePos,1);
    dim = size(coeffWave,2);
    
    T = 1/f;
    nt = floor(nCycles*T/dt);
    USources = zeros(nt,nSources,dim);
    t = (1:nt)'*dt;
    
    sigma = tau/10;
    envelop = 1./(1+exp(-(t-5*sigma)/sigma));
    

    %%
    absCoeff = abs(coeffWave)';
    phaseCoeff = angle(coeffWave)';
    
    for indDim = 1:dim
        USources(:,:,indDim) = sourceAmp*envelop.*sin(2*pi*f*t+phaseCoeff(indDim,:)).*absCoeff(indDim,:);
    end
    
    if 0 == 1
        %%
        newfig(1);
        subplot(1,2,1);
        imagesc(USources(:,:,1));
        colorbar;
        set(gca,'XTick',[]);
        set(gca,'YTick',[]);
        subplot(1,2,2);
        imagesc(USources(:,:,2));
        formatfig(1);
        set(gca,'XTick',[]);
        set(gca,'YTick',[]);
        colorbar;
    end
    
    paramSource.USources = USources;
    paramSource.sourcePos = sourcePos;
    paramSource.f = f;
    paramSource.nCycles = nCycles;
    paramSource.sourceAmp = sourceAmp;
    paramSource.forcing = forcing;
    
end










































