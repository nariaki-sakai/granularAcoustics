function [U0, paramLattice] = loadLattice(route,L,p,dk,strain,a,indFile)
    
    %[U0, paramLattice] = loadLattice(route,L,p,dk,strain,a,indFile)
    %       load sample lattice
    %
    %   INPUT
    %       route: root folder. Needs to be the folder that contains the
    %           subfolder Data and Codes
    %       L: sample lattice size, of the form L = (Ly,Lx)
    %       p: percolation parameter. Between 0 and 1
    %       dk: half-width of the uniform distribution centered at k=1 for
    %           random spring constants. dk=0 means no disorder from spring 
    %           strain: strain of the lattice
    %       a: length at rest of the springs
    %           indFile: file number of the sample. Usually it has 5 digits,
    %           the first two corresponds to the type of the sample (buffer size etc...) 
    %           and the last three are the indices of the sample among this type
    %
    %   OUTPUT
    %       U0: in plane displacement field at rest of the sample. Has size
    %           (Ly, Lx, 2) where the last dimension corresponds to
    %           displacement in x (1) or y (2)
    %       paramLattice: all parameters of the lattice and everything one
    %           would need is save in this structure class variable

    
    %%
    dirLattice = [route 'Data/lattices/'];
    
    fileName = sprintf('lattice_p%d%%dk%d%%L%dstrain%d%%springL%d%%',100*p,uint16(100*dk),L(1),100*strain,round(100*a));
    fileNb = num2str(indFile,'%5.0d');
    while length(fileNb) < 5
        fileNb = ['0' fileNb];
    end
    filePath = [dirLattice fileName '/' fileName '_' fileNb '.mat'];
    S = load(filePath);
    
    U0 = S.U0;
    paramLattice = S.paramLattice;
    
    
    %the different case below were for version compatibility during the
    %code developpement
    try
        paramSD = S.paramSD;
        detectorPos = paramSD.detectorPos;
        if isa(detectorPos,'double')
            tmpD = detectorPos;
            detectorPos = cell(2,1);
            detectorPos{1} = squeeze(tmpD(1,:,:));
            detectorPos{2} = squeeze(tmpD(2,:,:));
            
        end
        paramLattice.detectorPos = detectorPos;
        paramLattice.uniSource = paramSD.uniSource;
        save(filePath,'U0','paramLattice');
    catch
        detectorPos = paramLattice.detectorPos;
        if isa(detectorPos,'double')
            tmpD = detectorPos;
            detectorPos = cell(2,1);
            detectorPos{1} = squeeze(tmpD(1,:,:));
            detectorPos{2} = squeeze(tmpD(2,:,:));
            paramLattice.detectorPos = detectorPos;
            save(filePath,'U0','paramLattice');
        else
            paramLattice.detectorPos{1} = squeeze(paramLattice.detectorPos{1});
            paramLattice.detectorPos{2} = squeeze(paramLattice.detectorPos{2});
        end
    end
    try paramLattice.indSample;
    catch
        paramLattice.p = p;
        paramLattice.dk = dk;
        paramLattice.indSample = indFile;
        save(filePath,'U0','paramLattice');
    end
    
    if indFile ~= paramLattice.indSample
        paramLattice.indSample = indFile;
    end
    
    try paramLattice.a;
    catch
        paramLattice.a = a;
    end
    save(filePath,'U0','paramLattice');
    
        %Fix version compatibility
%     [paramLattice,paramSolver,paramSD] = fixVersionCompatibilities(paramLattice,paramSolver,paramSD);

end



















