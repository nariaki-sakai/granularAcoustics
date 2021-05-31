close all;
clearvars;

route = '../../../';
addpath(genpath([route 'Codes/']));

m0 = 1/(2*pi)^2;

%% PARAMETERS

%list of parameters for which we want to compute the scattering matrix. The
%list is a matix of size (nParam,2) where each row corresponds to a pair
%(p, dk) 
listParam = [0.85 0];           
nParam = size(listParam,1);

%LATTICE
L = 50*[1 1];               %size of the samples
strain = 0.2;               %strain
a = 1;                      %length of springs at rest

prefixSample = 21;          %tag of the sample type
tag = 'e-4';                %tag for naming the file for the scattering matrix.
                            %Here, 'e-4' means the scattering matrix has
                            %been computed with 

%SOURCE
f = 1.1;
optDisp = 0;

%SIMULATION
energyThreshold = 10^-4;
fRange = f+[-0.05 0.05];


nPool = 10;
try parpool(nPool)
catch
end

nLattices = 1;

dirData = sprintf('%sData/scatMatAllFreq/',route);

modeRange = guideModeIdxRange(f,L(1),strain);

%%
fprintf('\nPARAMETERS\n')
fprintf('f = %1.2f\n',f);
fprintf('Lx, Ly = [%d %d]\n',L(1),L(2));
fprintf('mode index min,max = [%1.2f %d]\n',modeRange(1),modeRange(2));


%% Create list of data to process


listSamples = zeros(nParam*nLattices,3);
for indParam = 1:nParam
    for indLattice = 1:nLattices
        listSamples((indLattice-1)*nParam+indParam,:) = [listParam(indParam,:) indLattice];
    end
end
nSamples = nParam*nLattices;


%% LOOP

for indParam = 1:nSamples
    
    p = listSamples(indParam,:);
    dk = p(2);
    indSample = p(3);
    p = p(1);
    
    fprintf('PARAM p=%1.2f dk=%1.2f sample n°%d\n',p,dk,indSample);
    
    dirScatMat = [dirData 'scatMatSH_' genDataName(L,p,dk,strain,a,prefixSample,f,'') '/'];
    createDir(dirScatMat);
    
    %%
    indFile = prefixSample*1000+indSample;
    charLattice = genDataName(L,p,dk,strain,a,indFile,f,tag);
    dataPath = sprintf('%sscatteringMatricesSH_%s.mat',dirScatMat,charLattice);
    
    %%
    if exist(dataPath,'file')
        fprintf('%s already exists\n\n',dataPath);
    else

        %%%%%%%%%
        %%LOAD
        %%%%%%%%%
        fprintf('SAMPLE %d / %d\n',indParam,size(listSamples,1));
        try
                [U0, paramLattice] = loadLattice(route,L,p,dk,strain,a,indFile);
        catch
            error('The lattice doesn''t exist');
        end


        %%
        %%%%%%%%%
        %COMPUTE
        %%%%%%%%%
        tic;
        [scatMat, fProbe] = scatteringMatrixPulseSHengine6(U0,paramLattice,f,energyThreshold,fRange,modeRange,optDisp);
        toc;

        %%
        %%%%%%%%%
        %SAVE
        %%%%%%%%%
        fprintf('save %s \n',dataPath);
        save(dataPath,'scatMat','fProbe');
        fprintf('done\n\n');
        
    end
    
end




disp('END')



















