clear all;
close all;

route = '../../..';
addpath(genpath(route));

dirLattices = [route '/Data/lattices/'];


%% PARAMETERS

%RANDOMNESS
p = 0.85;                   %percolation parameter
dk = 0;                     %spring constant disorder
strain = 0.2;               %strain
a = 1;                      %length at rest of springs
nNN = 4;                    %number of neighbours. Needs to be 4 or 8

%LATTICE PARAMETERS
L = 50*[1 1];               %size (Ly,Lx) of the sample
LAbsLayer = 250;            %length of the absorbing layer
Lbuffer = 50;               %length of the buffer
LtotLayer = LAbsLayer + Lbuffer;
nSamples = 1;               %Number of samples we want to create
prefixSample = '21';        %tag of the samples

%ALGORITHM PARAMETERS
%stop criterium: the motion of masses during 1000 time steps does not
%exceed dispThreshold (in lattice size unit)
dispThreshold = 10^-5;      
%the number of strain steps. At the limit of infinite strainStep, we have a
%quasistatic strain. dStrain is thus the strain that is applied at each
%stage of the relaxation
nStrain = 5;
dStrain = strain/nSTrain;
%Each intermediate step of the strain relaxation stops when the motion of
%all masses does not exceed quasistaticDispThre. Usually we don't need a
%value as low as dispThreshold 
quasistaticDispThre = 0.05;


%DISPLAY PARAM
fprintf('p = %1.2f - dk = %1.2f\n',p,dk);
disp(['Sample size = ' num2str(L) ]);
disp(['Size buffer = ' num2str(Lbuffer)]);
disp(['Size absorbing layer = ' num2str(LAbsLayer)]);



%% OUTPUT FILES

latticeName = sprintf('lattice_p%d%%dk%d%%L%dstrain%d%%springL%d%%',100*p,uint16(100*dk),L(1),100*strain,round(100*a));
subDir = [dirLattices latticeName '/'];
createDir(subDir);

%DATA TO SAVE
U0all = cell(nSamples,1);
paramLatticeAll = cell(nSamples,1);


%% LOOP
parfor indSample = 1:nSamples

    
    %% GENERATE RANDOM LATTICE

    fprintf('\nGENERATE RANDOM LATTICE\n');
    tic;
    [U0, paramLattice, EtotVSt] = randLattice(L,Lbuffer,LAbsLayer,p,dk,nNN,strain,a,dispThreshold,quasistaticDispThre,nStrain);
    toc;

    paramLatticeAll{indSample} = paramLattice;
    U0all{indSample} = U0;

end

%%
disp('SAVE');
for indSample = 1:nSamples
    
    U0 = U0all{indSample};
    paramLattice = paramLatticeAll{indSample};
    
    indFile = 1;
    fileId = num2str(indFile,'%03d');
    fileId = [num2str(prefixSample,'%02d') fileId];
    filePath = [subDir latticeName '_' fileId '.mat'];
    while exist(filePath,'file')
        indFile = indFile + 1;
        fileId = num2str(indFile,'%03d');
        fileId = [num2str(prefixSample,'%02d') fileId];
	filePath = [subDir latticeName '_' fileId '.mat'];
    end

    paramLattice.indSample = indFile;
    fprintf('file nb %d\n',indFile);
    save(filePath,'U0','paramLattice');
    
end

disp('END');













































