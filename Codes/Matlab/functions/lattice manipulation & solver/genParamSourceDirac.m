function paramSource = genParamSourceDirac(sourceAmp,coeffWave,sourcePos,forcing,lenDirac)

    %%

    nSources = size(sourcePos,1);
    dim = length(sourceAmp);
    
    USources = zeros(lenDirac,nSources,dim);
    for indDim = 1:dim
        USources(1,:,indDim) = sourceAmp(indDim)*coeffWave;
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
    paramSource.f = 'dirac';
    paramSource.sourceAmp = sourceAmp;
    paramSource.forcing = forcing;
    
    
end










































