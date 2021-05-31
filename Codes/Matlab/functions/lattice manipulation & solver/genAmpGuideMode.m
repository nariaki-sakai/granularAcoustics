function ampWave = genAmpGuideMode(modeIn,Ly,waveType)

    %ampWave = genAmpGuideMode(modeIn,L,waveType)
    %   generate the wave front amplitude of the source signal
    %
    %   INPUT:
    %       modeIn: either the index of the single mode that we want to
    %           send, or a vector of the complex amplitudes for a
    %           source signal that is a superposition of modes. In the latter
    %           case, the elements of the vector are the amplitudes for mode 1
    %           to mode length(modeIn)
    %       Ly: number of mass in the width the guide. The width of the
    %           guide thus equal Ly-1
    %           waveType: char which is either P or S, for compression of shear
    %           wave
    %
    %   OUTPUT:
    %       ampWave: vector corresponding to wave front amplitude of each
    %           mass of the source. ampWave has length equal to Ly

    z = linspace(0,Ly-1,Ly)';

    if length(modeIn) == 1
        switch waveType
            case 'P'
                ampWave = cos(modeIn*pi/(Ly-1)*z);
            case 'S'
                ampWave = sin(modeIn*pi/(Ly-1)*z);
            otherwise
                error('waveType has wrong attribute');
        end
    else
        modeIn  = reshape(modeIn,length(modeIn),1);
        modeIdx = (1:length(modeIn))';
        modeAmp = abs(modeIn);
        modePhase = angle(modeIn);
        
        switch waveType
            case 'P'
                ampWave = modeAmp.*cos(modeIdx*pi/(Ly-1)*z + modePhase);
            case 'S'
                ampWave = modeAmp.*sin(modeIdx*pi/(Ly-1)*z + modePhase);
            otherwise
                error('waveType has wrong attribute');
        end
        ampWave = sum(ampWave,1);
        
        
    end
    

    
end






















