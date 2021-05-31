close all;
clear all;

route = '../../../';
addpath(genpath(route));


%%

modeIdx = 4;

%lattice
L = 20*[1 1];
p = 0.85;
dk = 0;
strain = 0.2;
a = 1;


prefixSample = 31;
indSample = 1;
indFile = prefixSample*1000+indSample;

%signal sources and detection
sourceAmp = 0.1;

%solver parameters
periodSampling = 50;
freqSave = 1;
nIter = 1000;

frictionCoeff = 0.00;

forcing = 'position';
optBC = {'free';'clamped';'clamped';'free'};
waveType = 'S';

disp('LOAD LATTICE');
[U0, paramLattice] = loadLattice(route,L,p,dk,strain,a,indFile);
[U0, paramLattice] = isolateSample(U0, paramLattice, 0);
L = paramLattice.L;

paramLattice.detectorPos = [];
friction = paramLattice.friction;
paramLattice.friction = frictionCoeff/(2*pi)*ones(size(friction));

% f = 0;
% if f > 0
%     waveVec = waveVecGuideMode(f,paramLattice);
%     waveVec = waveVec(modeIdx);
% else
%     waveVec = 0;
% end
% 
% m = 1/(2*pi)^2;
% c = 2*pi*sqrt(strain/(1+strain));
% rho = 1/(2*pi)^2;
% Z0 = rho*c;
% Z = Z0*exp(-1i*waveVec/2);
% % Z = abs(Z);
% waveVec/2/pi;
% Z;

m = 1/(2*pi)^2;
c = 2*pi*sqrt(strain/(1+strain));
rho = 1/(2*pi)^2;
Z0 = rho*c;

Z = Z0;

imZ = -0;


viscousCoef = Z;
viscous = paramLattice.viscous;
viscous(:,1) = viscousCoef;
viscous(:,end) = viscousCoef;
paramLattice.viscous = viscous;

deltaMass = imZ*m;
newMass = m + deltaMass;
mMat = paramLattice.mMat;
mMat(:,1)   = newMass;
mMat(:,end) = newMass;
paramLattice.mMat = mMat;



%%


nIterFactor = 200;
%%%%%%%%%
%%LOAD
%%%%%%%%%
f = [];
sourcePos = [(1:L(1))' L(2)*ones(L(1),1)];
paramSolver = genParamSolver(periodSampling,nIterFactor,freqSave,f,paramLattice,optBC);
paramSolver.nIter = 1000000;
% paramSolver.optOut = 'sample';
paramSolver.optOut = 'detector';




lenDirac = 250;
coeffWave = genAmpGuideMode(modeIdx,L(1),'S');
coeffWave = zeros(size(coeffWave));
coeffWave(10) = 1;
paramSource = genParamSourceDirac(sourceAmp,coeffWave,sourcePos,forcing,lenDirac);


%%%DISPLAY PARAMETERS
displayParam(paramLattice,paramSource,paramSolver);
nDetector = 49;
detectorPos = [randi(L(1),nDetector,1) randi(L(2),nDetector,1)];
paramLattice.detectorPos = {detectorPos};


%% SOLVE EQUATIONS

disp('PROPAGATE WAVE');
tic;
Ui = zeros(size(U0));
Ui = Ui(:,:,1);
dUi = Ui;
dataOut = verletSolverSH(Ui,dUi,U0,paramLattice,paramSolver,paramSource);
toc;

UvsTime = dataOut.U;


%% ENERGY

time = dataOut.time;
EvsTime = dataOut.E;
newfig(1);
subplot(1,2,1);
hold all;
EvsTime(:,2) = EvsTime(:,2) - EvsTime(1,2);
plot(time,EvsTime(:,1),'r');
plot(time,EvsTime(:,2),'b');
plot(time,sum(EvsTime,2),'k');
xlabel('time');
ylabel('energy');
legend('E_c','E_p','E_{tot}','fontsize',12,'Location','NorthEast');
% set(gca,'YScale','log');
set(gca,'fontsize',12);

subplot(1,2,2);
hold all;
plot(time,EvsTime(:,1),'r');
plot(time,EvsTime(:,2),'b');
plot(time,sum(EvsTime,2),'k');
xlabel('time');
% ylabel('energy');
legend('E_c','E_p','E_{tot}','fontsize',12,'Location','SouthEast');
set(gca,'YScale','log');
% formatfig(1);
set(gca,'fontsize',12);


%% VIDEO
    
charLattice = genDataName(L,p,dk,strain,a,indSample,0,sprintf('realZ%1.3fimZ%1.3fm',Z,imZ));
if 0 == 1
    %%
    dirVideo = [route 'video/propagateDiracDoubleLeakComplex/'];
    createDir(dirVideo);


    vidPath = [dirVideo 'dispFieldZ_' charLattice];
    paramVideo = genParamVideo(paramSolver);
    paramVideo.freqFrame = 20 ;
    paramVideo.UScale = 50;
    paramVideo.tRange = [0 100];
    paramVideo.optOut = 'sample';

    dispVideo(2,dataOut,1,'position',paramLattice,paramVideo,paramSource);
end

%% VIBRATIONAL MODES

optBC = {'open','clamped','clamped','open'};
[eigenfreqs, eigenvecs, meshxD, meshyD] = vibrationalModes(U0,paramLattice,optBC,Z0);


%%
Lx = L(2)-1;
Ly = L(1)-1;
nx = 1:L(2)-2;
ny = (1:2:L(1)-2)';



kx = (nx-1/2)*pi/(Lx+1/2);
kx = nx*pi/Lx;
ky = ny*pi/Ly;
sinx = sin(kx/2);
siny = sin(ky/2);
s = paramLattice.strain;
f0 = (s/(1+s))^(1/2);

f_res =2*f0*sqrt(sinx.^2 + siny.^2);
[meshx,meshy] = meshgrid(nx,ny);
[f_res, ind] = sort(f_res(:));
meshx = meshx(ind);
meshy = meshy(ind);

eigenfreqs0 = f_res;




%% REAL SPAC

newfig(2);
imagesc(mMat);
axis equal;
xlim([0 L(2)]+1/2);
ylim([0 L(1)]+1/2);
for indD = 1:nDetector
    plot(detectorPos(indD,2),detectorPos(indD,1),'r.','MarkerSize',10);
end


%%

indDetector = 1;
signalOut = dataOut.U;
time = dataOut.time;
nt = length(time);
signalOut = signalOut(:,indDetector);

newfig(3);
plot(time,signalOut);


%% FOURIER SPACE

[fFFT, fftSignalOut] = fourierTransform(time,signalOut,nt,'energy');
% [fFFT, fftSignalIn ] = fourierTransform(time,signalIn ,nt,'energy');

newfig(4);
hold all;
plot(2*pi*fFFT,abs(fftSignalOut));
% plot(fFFT,abs(fftSignalIn));

% 
% set(gca,'YScale','log');
% set(gca,'XScale','log');
xlim([0 8]);
% xlim([0.1 2]);
% ylim([10^-5 0.1]);
% plot(real(eigenfreqs),10^-5 ,'r+');
% plot(real(eigenfreqs),abs(imag(eigenfreqs))./abs(real(eigenfreqs))*10^12,'r+');
formatfig([1000 200]);
xlabel('2\pi f');
ylabel(sprintf('fft[s_{x=%d,y=%d}]',detectorPos(indDetector,2),detectorPos(indDetector,1)));

h = figure(4);
winPos = h.Position;
winPos(1:2) = 0;
h.Position = winPos;

% strX = string(meshx);
% strY = string(meshy);
% str = "("+strX+","+strY+")";
% text(2*pi*eigenfreqs0-0.05,10^-4*ones(size(eigenfreqs0)),str);
if imag(Z) < 0
    title(sprintf('Z_0 = %1.3f%1.3fi',real(Z),imag(Z)));
else
    title(sprintf('Z_0 = %1.3f+%1.3fi',real(Z),imag(Z)));
end
grid on;



%%
newfig(5);
hold all;
plot(2*pi*fFFT,abs(fftSignalOut));
% plot(fFFT,abs(fftSignalIn));

% 
% set(gca,'YScale','log');
% set(gca,'XScale','log');
xlim([6.2 6.8]);
% xlim([0.1 2]);
% ylim([10^-5 0.1]);
% plot(real(eigenfreqs),0*10^-4 ,'r+');
% plot(real(eigenfreqs),abs(imag(eigenfreqs))./abs(real(eigenfreqs))*10^12,'r+');
formatfig([1000 200]);
xlabel('2\pi f');
ylabel(sprintf('fft[s_{x=%d,y=%d}]',detectorPos(indDetector,2),detectorPos(indDetector,1)));

h = figure(5);
winPos = h.Position;
winPos(1:2) = 0;
h.Position = winPos;

% strX = string(meshx);
% strY = string(meshy);
% str = "("+strX+","+strY+")";
% text(2*pi*eigenfreqs0-0.05,10^-4*ones(size(eigenfreqs0)),str);
if imag(Z) < 0
    title(sprintf('Z_0 = %1.3f%1.3fi',real(Z),imag(Z)));
else
    title(sprintf('Z_0 = %1.3f+%1.3fi',real(Z),imag(Z)));
end
grid on;

%%

newfig(6);
omega = 6.5;
ind = findClosest(omega,real(eigenfreqs));

vibModeField = zeros(L);
for indP = 1:length(meshxD)
    indi = meshyD(indP);
    indj = meshxD(indP);
    vibModeField(indi,indj) = abs(eigenvecs(indP,ind));
end
% imagesc(log10(vibModeField),[-5 0]);
imagesc(vibModeField);
xlim([0 L(2)]+0.5);
ylim([0 L(1)]+0.5);
title(sprintf('f = %1.2f + i%1.2e',real(eigenfreqs(ind)),imag(eigenfreqs(ind))));
axis equal;
colorbar

%% SAVEFIG

if 0 == 1
    %%
    
    dirFig = [route 'figs/propagateDiracClosedSystemLeakDoubleComplexZPointSource/' charLattice '/'];
    createDir(dirFig);
    
    figure(4);
    figPath = [dirFig 'fftSignal'];
    save2pdf(figPath,4);
    savefig(4,figPath);
    
    figure(5);
    figPath = [dirFig 'fftSignalZoom'];
    save2pdf(figPath,5);
    savefig(5,figPath);

    figure(6);
    figPath = [dirFig 'vibMode'];
    save2pdf(figPath,6);
    savefig(6,figPath);
    
end

disp('END');




















