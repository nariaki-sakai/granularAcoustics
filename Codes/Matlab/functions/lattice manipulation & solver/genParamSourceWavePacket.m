function paramSource = genParamSourceWavePacket(f,sourceAmp,nCycles,tau,U0,shapeWave,sourcePos,forcing,paramSolver)

    %paramSource = genParamSourceWavePacket(f,sourceAmp,nCycles,tau,U0,shapeWave,sourcePos,forcing,paramSolver)
    %   generates the structure class variable paramSource that contains
    %   all information one would need concerning the source signal. This
    %   function if for a wave packet; if one would need other type of
    %   source (Dirac, monochromatic), one should look at
    %   genParamSourceMonoChrome and genParamSourceDirac
    %
    %   INPUT:
    %       f: central frequency of the wave packet (gaussian)
    %       sourceAmp: amplitude of the gaussian of the wave packet
    %       nCycles: duration of the signal in number of period
    %       tau: size of the wave packet. The half-height width of the
    %           gaussian is tau/8, this gives a gaussian which magnitude is
    %           ~10^-3 at the beginning and end of the signal
    %       U0: in-plane dispalcement field of the lattice at rest
    %       shapeWave: shape of the wave front of the source signal, output
    %           by genAmpGuideMode
    %       sourcePos: position of the source. Is a (Ly,2) size matrix
    %           where each row is the indices (i,j) of the source
    %       forcing: string of character that is either 'force' or
    %           'position'; for the type of excitation
    %       paramSolver: the structure output by genParamSolver that
    %           contains all quantity one would need to solve numerically the
    %           dynamical equations
    
    
    if f > 0
        %%
        L = size(U0);
        dt = paramSolver.dt;

        nSources = size(sourcePos,1);
        dim = size(shapeWave,2);

        T = 1/f;
        nt = floor(nCycles*T/dt);
        USources = zeros(nt,nSources,dim);
        t = (1:nt)'*dt;

        sigma = tau/8;
        envelop = exp(-(t-tau/2).^2/(2*sigma^2));
        
        
        if 0 == 1
            %%
            newfig(13);
            plot(envelop);
        end

        %%
        absCoeff = abs(shapeWave)';
        phaseCoeff = angle(shapeWave)';
        
        globalAmp   = abs(sourceAmp);
        globalPhase = angle(sourceAmp);

        for indDim = 1:dim
            USources(:,:,indDim) = globalAmp*envelop.*sin(2*pi*f*t+phaseCoeff(indDim,:)+globalPhase).*absCoeff(indDim,:);
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
    
    else
        error('the frequency is zero');
    end
    
    paramSource.USources = USources;
    paramSource.sourcePos = sourcePos;
    paramSource.f = f;
    paramSource.nCycles = nCycles;
    paramSource.sourceAmp = sourceAmp;
    paramSource.forcing = forcing;
    paramSource.nCycles = nCycles;
    paramSource.tau = tau;
    
    
end










































