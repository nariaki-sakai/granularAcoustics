function waveVec = waveVecGuideMode(varargin)

    fSampling = varargin{1};
    paramLattice = varargin{2};
    switch nargin
        case 2
            nmin = 1;
            nmax = maxGuideModeIdx(max(fSampling),paramLattice);
        case 3
            modeRange = varargin{3};
            nmin = ceil(modeRange(1));
            nmax = floor(modeRange(2));
    end
    %%
    if isrow(fSampling)
        fSampling = fSampling';
    end
    
    strain = paramLattice.strain;
    omega0 = 2*pi*sqrt(strain/(1+strain));
    L = paramLattice.L;
    L = L(1)-1;
    
    n = nmin:nmax;
    omega = 2*pi*fSampling;
    waveVec = 2*asin(sqrt( (omega/(2*omega0)).^2 - sin(n*pi/(2*L)).^2 ));

end