close all;
clearvars;

%%
filePath = '/home/nariaki/Works/granularAcoustics/granularAcoustics_master/Data/scatMatAllFreq/scatMatSH_L20p85%dk0%strain20%springL100%sample31f1_050000e+00/scatteringMatricesSH_L20p85%dk0%strain20%springL100%sample31001f1_050000e+00e-4.mat';
load(filePath);



newfig(1);
nf = length(fProbe);
nModes = size(scatMat,1)/2;

S = scatMat(:,:,(nf-1)/2);
imagesc(abs(S*S'));

tf = scatMat(nModes+1:end,1:nModes,:);

Tf = zeros(nf,nModes);
for indf = 1:nf
    tmpTf = svd(tf(:,:,indf));
    Tf(indf,:) = tmpTf';
end

%%
newfig(2);
ylim([0 1]);
for indf = 1:10:nf
    cla;
    plot(Tf(indf,:));
    title(sprintf('f = %1.3f',fProbe(indf)));
    drawnow;
end


%%
newfig(3);

plot(2*pi*fProbe,Tf(:,1));
plot(2*pi*fProbe,Tf(:,2));
plot(2*pi*fProbe,Tf(:,3));
plot(2*pi*fProbe,sum(Tf,2));
legend('T_1','T_2','T_3','total')


