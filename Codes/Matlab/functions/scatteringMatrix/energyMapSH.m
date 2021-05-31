function [Ec, Ep] = energyMapSH(U,dU,U0,paramLattice)

    %%
    %%%INPUT
    kMat = paramLattice.kMat;
    mMat = paramLattice.mMat;
    m0 = 1/(2*pi)^2;
    a = paramLattice.a;
    LAbsLayer = paramLattice.LAbsLayer;
    LtotLayer = paramLattice.LtotLayer;
    Lbuffer = paramLattice.Lbuffer;
    L = paramLattice.L;

    distNode = distanceNode(U,U0);

    %%%ENERGY
    Ec = m0*sum(dU.^2,3)/2;
    Ep = 1/2*sum(kMat.*(distNode - a).^2,3)/2;
    Ec(mMat == 0) = 0;
    Ep(mMat == 0) = 0;

    Ec = Ec(:,LAbsLayer+1:LtotLayer+L(2)+Lbuffer);
    Ep = Ep(:,LAbsLayer+1:LtotLayer+L(2)+Lbuffer);
%     
%     kMat = paramLattice.kMat;
%     mMat = paramLattice.mMat;
%     m0 = 1/(2*pi)^2;
%     a = paramLattice.a;
% 
%     distNode = distanceNode(U,U0);
%     
%     %%%ENERGY
%     Ec = m0*dU.^2/2;
%     Ep = 1/2*sum(kMat.*(distNode - a).^2,3)/2;
%     Ec(mMat == 0) = 0;
%     Ep(mMat == 0) = 0;


end

function distNode = distanceNode(U,U0)

    Ux = U0(:,:,1);
    Uy = U0(:,:,2);
    Uz = U;
    Unx = neighbourVal(Ux);
    Uny = neighbourVal(Uy);
    Unz = neighbourVal(Uz);
    
    [distNodeX, distNodeY] = distanceNodeXY(Ux,Uy,Unx,Uny);
    distNodeZ = distanceNodeZ(Uz,Unz);
    
    distNode = sqrt(distNodeX.^2 + distNodeY.^2 + distNodeZ.^2);

end

function distNodeZ = distanceNodeZ(Uz,Unz)

    distNodeZ = zeros(size(Unz));
    
    vecLattice = [-1 -1;
               0 -1;
               1 -1;
              -1  0;
               1  0;
              -1  1;
               0  1;
               1  1];
           
    for indNN = 1:4
        distNodeZ(:,:,indNN) = Unz(:,:,indNN) - Uz;
    end
    for indNN = 5:8
        distNodeZ(:,:,indNN) = -circshift(distNodeZ(:,:,9-indNN),vecLattice(9-indNN,:));
    end
        
end

function [distNodeX, distNodeY] = distanceNodeXY(Ux,Uy,Unx,Uny)

    distNodeX = zeros(size(Unx));
    distNodeY = zeros(size(Uny));
    
    vecLattice = [-1 -1;
               0 -1;
               1 -1;
              -1  0;
               1  0;
              -1  1;
               0  1;
               1  1];
    for indNN = 1:4
        distNodeX(:,:,indNN) = Unx(:,:,indNN) + vecLattice(indNN,2) - Ux;
        distNodeY(:,:,indNN) = Uny(:,:,indNN) + vecLattice(indNN,1) - Uy;
    end
    for indNN = 5:8
        distNodeX(:,:,indNN) = -circshift(distNodeX(:,:,9-indNN),vecLattice(9-indNN,:));
        distNodeY(:,:,indNN) = -circshift(distNodeY(:,:,9-indNN),vecLattice(9-indNN,:));
    end
        
end

function Un = neighbourVal(U)

    %%
    vecLattice = [-1 -1;
           0 -1;
           1 -1;
          -1  0;
           1  0;
          -1  1;
           0  1;
           1  1];

    %%% X
    Un = zeros([size(U) 8]);
    for ind = 1:8
        Un(:,:,ind) = circshift(U,-vecLattice(ind,:));
    end

end




















