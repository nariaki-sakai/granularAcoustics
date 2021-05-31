function [eigenfreqs, eigenvecs, meshxD, meshyD] = vibrationalModes(varargin)

    m = 1/(2*pi)^2;

    U0 = varargin{1};
    paramLattice = varargin{2};
    optBC = varargin{3};
    
    switch nargin
        case 3
            Z = 0;
        case 4
            Z = varargin{4};
    end
    
    %%%
    [D, V, meshxD, meshyD] = dynamicMatrixSH(U0,paramLattice,optBC);
    
    tic;
    disp('DIAGONALIZATION: might take a while...');
    [eigenvecs,eigenvals]=polyeig(D,Z/m*V,-eye(size(D)));
    disp('DONE');
    toc;
    eigenfreqs = imag(eigenvals) + 1i*real(eigenvals);
    

end


function [D, V, meshx, meshy] = dynamicMatrixSH(U0,paramLattice,optBC)

    %assemble_dynamic_matrix returns the Dynamic matrix D in order to compute
    %the normal modes of the network
    %   D is the block matrix : (A,Gamma;Gamma,B)
    %   K_mat is a n*n matrix and stores the stiffnesses between mass i and mass j
    %(symmetric)
    %   M_mat in a n*n matrix and stores the masses mij=mi. M_mat is in general
    %non-symmetric, except if all the masses are equal to the same value.
    %   L0_mat is a n*n matrix and stores the length at rest of spring between
    %   masses i and j. This is a symmetric matrix.
    %   L_eq_x_mat is a n*n matrix and stores the x components of the distance
    %   vector going from mass i to mass j. This is an antisymmetric matrix
    %   L_eq_y_mat is a n*n matrix and stores the y components of the distance
    %   vector going from mass i to mass j. This is an antisymmetric matrix
    
    mMat = paramLattice.mMat;
    kMat0 = paramLattice.kMat;
%     strain = paramLattice.strain;
    a = paramLattice.a;
    L = size(mMat);
        
    if ~isempty(kMatCheck(kMat0))
        error('kMat not symmetric 1');
    end

    kMat = zeros(prod(L));
    for indP = 1:prod(L)
        shiftVect = [-1 -1;
                      0 -1;
                      1 -1;
                     -1 0;
                      1 0;
                     -1 1;
                      0 1;
                      1 1]; 
                  
        for indShift = 1:8
            indN = neighbourId(indP,L,shiftVect(indShift,:));
            if ~isempty(indN)
                [indi,indj] = ind2sub(L,indP);
                kMat(indP,indN) = kMat0(indi,indj,indShift);
            end
        end
    end
    
    if 0 == 1
        for indP = 1:prod(L)
            if any(kMat(indP,:) ~= kMat(:,indP)')
                error('kMat not symmetric 2');
            end
        end
    end
    

    

    %% Boundary conditions and missing masses
    
    [meshx,meshy] = meshgrid(1:L(2),1:L(1));
    boundaryPos = [1 1 L(1) L(2)];
    meshXY = {meshx, meshy, meshy, meshx};
        %% 
    V = zeros(prod(L));
    for indBC = 1:4
        switch optBC{indBC}
            case 'clamped'
                mMat(meshXY{indBC} == boundaryPos(indBC)) = inf;
            case 'open'
                switch indBC
                    case 1
                        ind = find(meshXY{indBC}(:) == 1);
                    case 2
                        ind = find(meshXY{indBC}(:) == L(2));
                    case 3
                        ind = find(meshXY{indBC}(:) == 1);
                    case 4
                        ind = find(meshXY{indBC}(:) == L(1));
                end
                for ii = 1:length(ind)
                    V(ind(ii),ind(ii)) = 1;
                end
        end
    end
    masses_to_remove = find(mMat(:) == 0 | mMat(:) == inf);
    masses_to_remove = flipud(masses_to_remove);
    
        
    %% relative position between nodes n and p
    
    relPos = relativePos(U0);
    Leq   = sqrt(relPos(:,:,1).^2+relPos(:,:,2).^2);
    
    M_mat = ones(1,length(mMat(:)));
    M_mat = mMat(:)*M_mat;
    
    
    %% Dynamical matrix
    
    D = kMat ./ M_mat .*( 1 - a./Leq );
    D(isnan(D))=0;
    sum_cols_A=sum(D,2);
    D=D-diag(sum_cols_A);
    
    meshx = meshx(:);
    meshy = meshy(:);

    
    %%
    
    D(masses_to_remove,:) = [];
    D(:,masses_to_remove) = [];
    V(masses_to_remove,:) = [];
    V(:,masses_to_remove) = [];
    meshx(masses_to_remove) = [];
    meshy(masses_to_remove) = [];



    
end


function Leq = relativePos(U0)
    %L0_Leq_Mmat_Kmat_create Creates the matrices L0, Leq, Leq_x,Leq_y, given the
    %equilibrium state and length a rest a0 of the springs

    Ux0 = U0(:,:,1);
    Uy0 = U0(:,:,2);
    
    L = size(Ux0);
    
    [ lengths_down_left_x,lengths_down_left_y ,lengths_down_x,lengths_down_y,lengths_down_right_x,lengths_down_right_y, lengths_right_x,lengths_right_y,lengths_up_right_x,lengths_up_right_y,lengths_up_x,lengths_up_y,lengths_up_left_x,lengths_up_left_y,lengths_left_x,lengths_left_y,lengths_down_left,lengths_down,lengths_down_right,lengths_right,lengths_up_right,lengths_up,lengths_up_left,lengths_left] = lengths_vectors( Ux0,Uy0,1);

    Leq_x=zeros(prod(L));
    Leq_y=zeros(prod(L));

    for p=1:prod(L)

        [ row,col ] = ind2sub(L,p);
        if 1
            
            %% Finding neighbours
            neigh_dl=neighbourId(p,L,[-1 -1]);
            neigh_d=neighbourId(p,L,[-1 0]);
            neigh_dr=neighbourId(p,L,[-1 1]);
            neigh_l=neighbourId(p,L,[0 -1]);
            neigh_r=neighbourId(p,L,[0 1]);
            neigh_ul=neighbourId(p,L,[1 -1]);
            neigh_u=neighbourId(p,L,[1 0]);
            neigh_ur=neighbourId(p,L,[1 1]);



            %% Leq_x
            Leq_x(p,neigh_dl)=lengths_down_left_x(row,col ) ;
            Leq_x(neigh_dl,p)=-Leq_x(p,neigh_dl);

            Leq_x(p,neigh_d)=lengths_down_x(row,col );
            Leq_x(neigh_d,p)=-Leq_x(p,neigh_d);

            Leq_x(p,neigh_dr)=lengths_down_right_x(row,col );
            Leq_x(neigh_dr,p)=-Leq_x(p,neigh_dr);

            Leq_x(p,neigh_r)=lengths_right_x(row,col );
            Leq_x(neigh_r,p)=-Leq_x(p,neigh_r);

            Leq_x(p,neigh_ur)=lengths_up_right_x(row,col );
            Leq_x(neigh_ur,p)=-Leq_x(p,neigh_ur);

            Leq_x(p,neigh_u)=lengths_up_x(row,col );
            Leq_x(neigh_u,p)=-Leq_x(p,neigh_u);

            Leq_x(p,neigh_ul)=lengths_up_left_x(row,col );
            Leq_x(neigh_ul,p)=-Leq_x(p,neigh_ul);

            Leq_x(p,neigh_l)=lengths_left_x(row,col );
            Leq_x(neigh_l,p)=-Leq_x(p,neigh_l);

            %% Leq_y
            Leq_y(p,neigh_dl)=lengths_down_left_y(row,col ) ;
            Leq_y(neigh_dl,p)=-Leq_y(p,neigh_dl);

            Leq_y(p,neigh_d)=lengths_down_y(row,col );
            Leq_y(neigh_d,p)=-Leq_y(p,neigh_d);

            Leq_y(p,neigh_dr)=lengths_down_right_y(row,col );
            Leq_y(neigh_dr,p)=-Leq_y(p,neigh_dr);

            Leq_y(p,neigh_r)=lengths_right_y(row,col );
            Leq_y(neigh_r,p)=-Leq_y(p,neigh_r);

            Leq_y(p,neigh_ur)=lengths_up_right_y(row,col );
            Leq_y(neigh_ur,p)=-Leq_y(p,neigh_ur);

            Leq_y(p,neigh_u)=lengths_up_y(row,col );
            Leq_y(neigh_u,p)=-Leq_y(p,neigh_u);

            Leq_y(p,neigh_ul)=lengths_up_left_y(row,col );
            Leq_y(neigh_ul,p)=-Leq_y(p,neigh_ul);

            Leq_y(p,neigh_l)=lengths_left_y(row,col );
            Leq_y(neigh_l,p)=-Leq_y(p,neigh_l); 

        end

    end

    Leq = zeros([prod(L) prod(L) 2]);
    Leq(:,:,1) = Leq_x;
    Leq(:,:,2) = Leq_y;

end

function notSymIdx = kMatCheck(kMat)

    notSymIdx = [];
    L = size(kMat(:,:,1));
    for ii = 2:L(1)-1
        for jj = 2:L(2)
            
            if ~(kMat(ii,jj,1) == kMat(ii-1,jj-1,8))
                notSymIdx = [notSymIdx; ii jj 1];
            end
            if ~(kMat(ii,jj,2) == kMat(ii,jj-1,7))
                notSymIdx = [notSymIdx; ii jj 2];
            end
            if ii < L(1)
                if ~(kMat(ii,jj,3) == kMat(ii+1,jj-1,6))
                    notSymIdx = [notSymIdx; ii jj 3];
                end
            end
            if ~(kMat(ii,jj,4) == kMat(ii-1,jj,5))
                notSymIdx = [notSymIdx; ii jj 4];
            end
        end
    end

    
end


function id = neighbourId(p,L,shiftVect)

    %%
    
    [indi,indj] = ind2sub(L,p);
    indi = indi + shiftVect(1);
    indj = indj + shiftVect(2);
    
    if indi > 0 & indi <= L(1) & indj > 0 & indj <= L(2)
        id = sub2ind(L,indi,indj);
    else
        id = [];
    end

end

function [ lengths_down_left_x,lengths_down_left_y ,lengths_down_x,lengths_down_y,lengths_down_right_x,lengths_down_right_y, lengths_right_x,lengths_right_y,lengths_up_right_x,lengths_up_right_y,lengths_up_x,lengths_up_y,lengths_up_left_x,lengths_up_left_y,lengths_left_x,lengths_left_y,lengths_down_left,lengths_down,lengths_down_right,lengths_right,lengths_up_right,lengths_up,lengths_up_left,lengths_left] = lengths_vectors( U_x,U_y,a )
%lengths_vectors Returns the vectors giving the lengths vectors (components
%and norms)

    [lengths_down_x,lengths_down_y]   = length_down(U_x,U_y,a);
    [lengths_right_x,lengths_right_y] = length_right(U_x,U_y,a);
    [lengths_down_right_x,lengths_down_right_y ] = length_down_right( U_x,U_y,a );
    [lengths_down_left_x,lengths_down_left_y ] = length_down_left( U_x,U_y,a );
    
    lengths_left_x = -circshift(lengths_right_x,[0,1]);
    lengths_left_y = -circshift(lengths_right_y,[0,1]);
    lengths_up_x   = -circshift(lengths_down_x,[-1,0]);
    lengths_up_y   = -circshift(lengths_down_y,[-1,0]);
    lengths_up_right_x = -circshift(lengths_down_left_x,[-1,-1]);
    lengths_up_right_y = -circshift(lengths_down_left_y,[-1,-1]);
    lengths_up_left_x  = -circshift(lengths_down_right_x,[-1,1]);
    lengths_up_left_y  = -circshift(lengths_down_right_y,[-1,1]);
    
    lengths_down  = sqrt(lengths_down_x.^2+lengths_down_y.^2);
    lengths_right = sqrt(lengths_right_x.^2+lengths_right_y.^2);
    lengths_up    = circshift(lengths_down,[-1,0]);
    lengths_left  = circshift(lengths_right,[0,1]);
    
    lengths_down_right=sqrt(lengths_down_right_x.^2+lengths_down_right_y.^2);
    lengths_down_left=sqrt(lengths_down_left_x.^2+lengths_down_left_y.^2);
    lengths_up_right=circshift(lengths_down_left,[-1,-1]);
    lengths_up_left=circshift(lengths_down_right,[-1,1]); 
    
    
    if any(isinf(lengths_up_left))
        error('lUL is inf');
    elseif any(isinf(lengths_left))
        error('lL is inf');
    elseif any(isinf(lengths_down_left))
        error('lDL is inf');
    elseif any(isinf(lengths_up))
        error('lU is inf');
    elseif any(isinf(lengths_down))
        error('lD is inf');
    elseif any(isinf(lengths_up_right))
        error('lUR is inf');
    elseif any(isinf(lengths_right))
        error('lR is inf');
    elseif any(isinf(lengths_down_right))
        error('lDR is inf');
    end
    
end

function [ l_down_x,l_down_y ] = length_down( Ux,Uy,a )
%length_down Returns a matrix with the lengths of down springs

l_down_x=circshift(Ux,[1,0])-Ux;
l_down_y=circshift(Uy,[1,0])-Uy-a;

end

function [ l_down_left_x,l_down_left_y ] = length_down_left( Ux,Uy,a )
%length_down_left Returns a matrix with the lengths of down-right (SE) springs

l_down_left_x=circshift(Ux,[1,1])-Ux-a;
l_down_left_y=circshift(Uy,[1,1])-Uy-a;

% figure
% subplot(1,2,1)
% imagesc(l_down_left_x)
% 
% subplot(1,2,2)
% imagesc(l_down_left_y)


end

function [ l_down_right_x,l_down_right_y ] = length_down_right( Ux,Uy,a )
%length_down_right Returns a matrix with the lengths of down-right (SE) springs

l_down_right_x=circshift(Ux,[1,-1])-Ux+a;
l_down_right_y=circshift(Uy,[1,-1])-Uy-a;

end

function [ l_right_x,l_right_y ] = length_right( Ux,Uy,a )
%length_right Returns a matrix with the lengths of right springs

l_right_x=circshift(Ux,[0,-1])-Ux+a;
l_right_y=circshift(Uy,[0,-1])-Uy;

end























































