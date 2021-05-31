function paramVideo = genParamVideo(paramSolver,varargin)
    
%%
    nIter = paramSolver.nIter;
    dt = paramSolver.dt;
    tRange = [0 nIter*dt];
    
    freqFrame = 1;
    zoomScale = 0;
    rattlers = 1;
    UScale = 1;
    
    paramVideo.tRange = tRange;
    paramVideo.freqFrame = freqFrame;
    paramVideo.zoomScale = zoomScale;
    paramVideo.rattlers = rattlers;
    paramVideo.UScale = UScale;
    paramVideo.optOut = paramSolver.optOut;
%     paramVideo.format = 'Motion JPEG AVI';
    paramVideo.format = 'MPEG-4';
    paramVideo.windowHeight = 600;
    paramVideo.quality = 100;
    paramVideo.frameRate = 20;
    paramVideo.cScale = 'lin';
    
    switch nargin
        case 1
            paramVideo.vidPath = [];
        case 2
            paramVideo.vidPath = varargin{1};
        otherwise
            error('number of inputs has to be 1 or 2');
    end

    
end
