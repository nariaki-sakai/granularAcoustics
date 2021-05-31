function [UOut, paramLatticeOut] = isolateSample(U0,paramLattice,expandSize)
    
    %%
    LtotLayer = paramLattice.LtotLayer;
    L = paramLattice.L;
    LAbsLayer = paramLattice.LAbsLayer;
    
    paramLatticeOut = paramLattice;
    
    mMat = paramLatticeOut.mMat;
    kMat = paramLatticeOut.kMat;
    viscous = paramLatticeOut.viscous;
    friction = paramLatticeOut.friction;
    
    UOut = U0(:,LtotLayer+1-expandSize:LtotLayer+L(2)+expandSize,:);
    mMat = mMat(:,LtotLayer+1-expandSize:LtotLayer+L(2)+expandSize);
    kMat = kMat(:,LtotLayer+1-expandSize:LtotLayer+L(2)+expandSize,:);
    viscous = viscous(:,LtotLayer+1-expandSize:LtotLayer+L(2)+expandSize);
    friction = friction(:,LtotLayer+1-expandSize:LtotLayer+L(2)+expandSize);

%     if expandSize
%         UOut = U0(:,LtotLayer+1-L(1):LtotLayer+L(2)+L(1),:);
%         mMat = mMat(:,LtotLayer+1-L(1):LtotLayer+L(2)+L(1));
%         kMat = kMat(:,LtotLayer+1-L(1):LtotLayer+L(2)+L(1),:);
%         viscous = viscous(:,LtotLayer+1-L(1):LtotLayer+L(2)+L(1));
%         friction = friction(:,LtotLayer+1-L(1):LtotLayer+L(2)+L(1));
%     else
%         UOut = U0(:,LtotLayer+1:LtotLayer+L(2),:);
%         mMat = mMat(:,LtotLayer+1:LtotLayer+L(2));
%         kMat = kMat(:,LtotLayer+1:LtotLayer+L(2),:);
%         viscous = viscous(:,LtotLayer+1:LtotLayer+L(2));
%         friction = friction(:,LtotLayer+1:LtotLayer+L(2));
%     end
    
    kMat(:,1,1:3) = 0;
    kMat(:,end,6:8) = 0;
    
    paramLatticeOut.mMat = mMat;
    paramLatticeOut.kMat = kMat;
    paramLatticeOut.viscous = viscous;
    paramLatticeOut.friction = friction;
    paramLatticeOut.LAbsLayer = 0;
    paramLatticeOut.LtotLayer = 0;
    paramLatticeOut.Lbuffer = 0;
    paramLatticeOut.detectorPos = [];
    
    if 0 == 1
        figure(13);
        subplot(1,2,1);
        cla;
        imagesc(mMat);
        subplot(1,2,2);
        cla;
        imagesc(kMat(:,:,2));
    end

end























