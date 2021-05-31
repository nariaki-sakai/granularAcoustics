function dispVideo(fig,data,indDir,optData,paramLattice,paramVideo,paramSource)

    %%

    
    
    switch optData
        case 'position'
            U = data.U;
            try
                U0 = data.U0;
                U = U - U0;
            catch
            end
        case 'velocity'
            U = data.dU;
        case 'flux'
            U = data.flux;
            U = U(:,:,:,1);
        otherwise
            error('optData has wrong attribute');
    end
    time = data.time;
    maxU = max(U(:));
    if maxU == 0
        warning('dataOut is zero');
        maxU = 1;
    end
    
    L = paramLattice.L;
    Lbuffer = paramLattice.Lbuffer;
    LAbsLayer = paramLattice.LAbsLayer;
    LtotLayer = paramLattice.LtotLayer;
    mMat = paramLattice.mMat;
    
    tRange = paramVideo.tRange;
    freqFrame = paramVideo.freqFrame;
    zoomScale = paramVideo.zoomScale;
    vidPath = paramVideo.vidPath;
    UScale = paramVideo.UScale;
    videoFormat = paramVideo.format;
    windowHeight = paramVideo.windowHeight;
    vidQuality = paramVideo.quality;
    cScale = paramVideo.cScale;
    
    detectorPos = paramLattice.detectorPos;
    sourceAmp = paramSource.sourceAmp;

    Ltot = size(U);
    Ltot = Ltot(2:3);
    
    try 
        leftPos  = detectorPos{1};
        rightPos = detectorPos{2};
    catch
        leftPos = [0 0];
        rightPos = [0 0];
    end
    leftPos  = squeeze(leftPos);
    rightPos = squeeze(rightPos);
    [iMass, jMass] = find(mMat == 0);

    screenSize = [0 0 1280 800];
    zoomScaleMax = 1 - 2*Ltot(1)/Ltot(2);
    if zoomScale > zoomScaleMax
        zoomScale = zoomScaleMax;
    end

    switch paramVideo.optOut
        case 'full'

        case 'part'
            if 2*Lbuffer+3*L(2) < L(2) + LAbsLayer + 2*Lbuffer
                leftPos(:,2)  = L(2)+1;
                rightPos(:,2) = Ltot(2) - L(2);
                jMass = jMass - (LAbsLayer - L(2));
            else
                leftPos(:,2)  = LAbsLayer/2+1;
                rightPos(:,2) = Ltot(2) - LAbsLayer/2;
                jMass = jMass - LAbsLayer/2;
            end
        case 'guide'
            try
                switch paramLattice.removeLead
                    case 'left'
                        leftPos(:,2)  = 1;
                        rightPos(:,2) = Ltot(2);
                        jMass = jMass - LtotLayer;
                    case 'right'
                        leftPos(:,2)  = 1;
                        rightPos(:,2) = Ltot(2);
                        jMass = jMass - LAbsLayer;
                end
            catch
                leftPos(:,2)  = 1;
                rightPos(:,2) = Ltot(2);
                jMass = jMass - LAbsLayer;
            end

        case 'sample'
            try
                switch paramLattice.removeLead
                    case 'left'
                    case 'right'
                        jMass = jMass - LtotLayer;
                end
            catch
                jMass = jMass - LtotLayer;
            end
            leftPos  = [];
            rightPos = [];
            screenSize = [0 0 400 400];
            zoomScale = paramVideo.zoomScale;
        case 'detector'
            leftPos  = [];
            rightPos = [];
        otherwise
            error('optOut has wrong attribute');
    end
    if isempty(detectorPos)
        leftPos  = [];
        rightPos = [];
    end
    

    %%
    if ~isempty(vidPath)
        if strcmp(computer,'GLNXA64') && strcmp(paramVideo.format,'MPEG-4')
            videoFormat = 'Motion JPEG AVI';
        end
        v = VideoWriter(vidPath,videoFormat);
        v.Quality = vidQuality;
        v.FrameRate = paramVideo.frameRate;
        open(v)
    end
    
    %%
    newfig(fig);
    
    axisPos = screenSize;
    axisPos(4) = ceil(axisPos(3)*Ltot(1)/Ltot(2)/(1-zoomScale));

    windowPos  = screenSize;
    windowPos(3) = axisPos(3);
    windowPos(1) = screenSize(3)/2-windowPos(3)/2;
    windowPos(2) = screenSize(4)/2-windowPos(4)/2;
    windowPos = windowPos*windowHeight/800;
    axisPos   = axisPos*windowHeight/800;
    windowPos(4) = axisPos(4) + 50;
    windowPos = floor(windowPos/2)*2;




    limX = Ltot(2)/2 + axisPos(3)/axisPos(4)*Ltot(1)/2*[-1 1];
    limX = round(limX);
    limX(1) = max(limX(1),1);
    limX(2) = min(limX(2),Ltot(2));
    
    mMat = mMat(:,limX(1):limX(2));


    [~, indt0] = min(abs(time - tRange(1)));
    if tRange(2) < inf
        [~, indt1] = min(abs(time - tRange(2)));
    else
        indt1 = length(time);
    end
    %%

    set(fig,'Units','pixels');
    set(fig,'Position',windowPos);
    set(fig,'Color',[1 1 1]);
    set(gca,'Units','pixels');
    set(gca,'Position',axisPos);
    
    %%
    for indt = indt0:freqFrame:indt1
        %%
        cla;
        hold all;
        tmpU = squeeze(U(indt,:,limX(1):limX(2),indDir));

        
        if sourceAmp > 0
            imagesc(tmpU,maxU/UScale*[-1 1]);
        else
            imagesc(tmpU,0.01*[-1 1]);
        end

        if ~(isempty(leftPos) || isempty(rightPos))
            if leftPos(1,2) >= limX(1) && rightPos(1,2) <= limX(2)
                plot(leftPos(:,2)-limX(1)+1,leftPos(:,1),'r.');
                plot(rightPos(:,2)-limX(1)+1,rightPos(:,1),'r.');
            end
        end
        figure(fig);
        plot(jMass-limX(1)+1,iMass,'k.');

        title(['t = ' num2str(time(indt),'%3.0f')]);
        axis equal;

        set(gca,'Units','pixels');
        set(gca,'Position',axisPos);
        set(fig,'Position',windowPos);
        drawnow;

        if ~isempty(vidPath)
            frame = getframe(fig);
            writeVideo(v,frame)
        end
    end

    if ~isempty(vidPath)
        close(v);
        
        %% CONVERT
        if strcmp(computer,'GLNXA64') && strcmp(paramVideo.format,'MPEG-4')
            
            aviPath = [vidPath '.avi'];
            mp4Path = [vidPath '.mp4'];
            [~,cmdout]=system(sprintf('ffmpeg -i %s -y -an -c:v libx264 -crf 20 -preset slow %s',aviPath,mp4Path));
            disp(cmdout);
            system(sprintf('rm %s',strrep(aviPath,'%','%')));
            
        end

    end
    

end











