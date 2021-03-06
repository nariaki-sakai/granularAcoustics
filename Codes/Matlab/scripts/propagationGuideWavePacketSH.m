close all;
clear all;

route = '../../../';
addpath(genpath(route));


%% PARAMETERS

%LATTICE PARAMETERS
L = 20*[1 1];           %Sample lattice size 
p = 0.85;               %percolation parameter
dk = 0;                 %parameter for random spring constant
strain = 0.2;           %strain
a = 1;                  %length at rest of the springs

indFile = 31001;        %fileNb

%SOURCE SIGNAL PARAMETERS
f = 1.04;               %excitation frequency
nCycles = 50;           %number of cycles in the wave packet
sourceAmp= 0.001;       %amplitude of the wave
tau = nCycles/f;        %duration of the source
modeIdx = 17;           %guide mode index of the excitation

%NUMERCIAL SOLVER PARAMETERS

%Number of sampling in one period of signal. Determines the time step dt
periodSampling = 50;    

%Determines the duration of the simulations, given by deltaT =
%nIterFactor*Ly/c, where Ly is the length of the wave guide, and c is a
%wave velocity chosen to be equal to c = 2*pi (which is actually an upper bound)
nIterFactor = 2;       

freqSave = 1;           %Frequency at which we save the state of the lattice.
waveType = 'S';         %Determines whether excitation is a shear wave or a compression wave.
forcing = 'position';   %Determines whether the excitation is made by imposing the position or the force.
sourceId = 1;           %Excitation from left (=1) or right (=2) of the guide.
optBC = 'clamped';      %Boundary conditions (BC) during the simulation. If one wish different BCs for different
                        %walls, optBC needs to be a (4,1) size cell where
                        %each element is a string giving the BC of each
                        %wall, beginning from left wall, bottom, top,
                        %right.



%%%%%%%%%%%%%%%%%%%%
%%PREPARE SIMULATION
%%%%%%%%%%%%%%%%%%%%

%LOAD LATTICE
[U0, paramLattice] = loadLattice(route,L,p,dk,strain,a,indFile);

%BUILD SOURCE SIGNAL AND SOLVER PARAMETER
sourcePos = paramLattice.detectorPos{sourceId};
ampWave = genAmpGuideMode(modeIdx,L(1),waveType);
paramSolver = genParamSolver(periodSampling,nIterFactor,freqSave,f,paramLattice,optBC);
paramSource = genParamSourceWavePacket(f,sourceAmp,nCycles,tau,U0,ampWave,sourcePos,forcing,paramSolver);
    
%CREATE TAG OF THE SIMULATION FOR SAVE
charLattice = genDataName(L,p,dk,strain,a,indFile,f,forcing);
charLattice = sprintf('%sMode%d',charLattice,modeIdx);

%DISPLAY PARAMETERS
displayParam(paramLattice,paramSource,paramSolver);





%% WAVE PACKET PROPAGATION

%Initial condition
Ui = zeros(size(U0,1),size(U0,2));  %position
dUi = Ui;                           %velocity

disp('PROPAGATE WAVE');
tic;
dataOut = verletSolverSH(Ui,dUi,U0,paramLattice,paramSolver,paramSource);
toc;



%% ENERGY

time = dataOut.time;
EvsTime = dataOut.E;
newfig(1);
EvsTime = EvsTime - EvsTime(1,:);
plot(time,EvsTime(:,1),'r');
plot(time,EvsTime(:,2),'b');
plot(time,sum(EvsTime,2),'k');
xlabel('time');
ylabel('energy');
legend('E_c','E_p','E_{tot}','fontsize',20,'Location','NorthEast');


%% DISPLAY

if 1 == 1
    %%
    
    dirVideo = sprintf('%svideo/dispFieldPulsesSH/L%dp%d%%dk%d%%strain%d%%restSL%d%%sample%d/',route,L(1),100*p,100*dk,100*strain,round(100*a),indFile);
    createDir(dirVideo);

    vidPath = [dirVideo 'dispFieldZ_' charLattice];
    paramVideo = genParamVideo(paramSolver);
    paramVideo.freqFrame = 100;
    paramVideo.UScale = 3;
    paramVideo.tRange = [0 inf];
    paramVideo.frameRate = 50;
    dispVideo(3,dataOut,1,'position',paramLattice,paramVideo,paramSource);
    
end


if 0 == 1
    %%
    
    vidPath = [dirVideo 'dispFieldZzoom_' charLattice];
    paramVideo = genParamVideo(paramSolver);
    paramVideo.freqFrame = 2;
    paramVideo.UScale = 1;
    paramVideo.zoomScale = 0.9;
    paramVideo.tRange = [60 160];
    dispVideo(4,dataOut,1,'position',paramLattice,paramVideo,paramSource);
    
end



if 0 == 1
    
    %%
    dirFig = [route 'figs/wavePropagationPulseSH/'];
    if ~exist(dirFig,'dir')
        mkdir(dirFig);
    end

    figure(1);
    figPath = [dirFig 'energy' charLattice];
    savefig(1,figPath);
    save2pdf(figPath,1)
    
    figure(3);
    figPath = [dirFig 'fftUmapIn' charLattice];
    savefig(3,figPath);
    save2pdf(figPath,3)
    
    figure(4);
    figPath = [dirFig 'fftUmapOut' charLattice];
    savefig(4,figPath);
    save2pdf(figPath,4)
    
    figure(5);
    figPath = [dirFig 'fftUlOut' charLattice];
    savefig(5,figPath);
    save2pdf(figPath,5)
    
    figure(6);
    figPath = [dirFig 'fftUrOut' charLattice];
    savefig(6,figPath);
    save2pdf(figPath,6)
    
    figure(7);
    figPath = [dirFig 'fftUlIn' charLattice];
    savefig(7,figPath);
    save2pdf(figPath,7)
    
    figure(8);
    figPath = [dirFig 'fftUrIn' charLattice];
    savefig(8,figPath);
    save2pdf(figPath,8)
    
end
disp('END')



















