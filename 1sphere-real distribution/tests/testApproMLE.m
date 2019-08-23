function [  ] = testApproMLE(  )

miu = 0;
sigma = 0.1;
theta0 = 0;
k = 1;
rho = [sqrt(0.95/2),sqrt(0.95/2)];
Ns = 100000;

NGrid = 100;

%% vary k
kV = logspace(-2,2,NGrid)';
PV = rho*sigma./sqrt(kV);

kDiff = zeros(NGrid,1);
theta0Diff = zeros(NGrid,1);
KLDiv = zeros(NGrid,1);
parfor ng = 1:NGrid
    [x,theta] = cylinderSampling(miu,sigma,PV(ng,:),kV(ng),theta0,Ns);
    [miuE,sigmaE,PE,kE,theta0E] = cylinderDistOneParaMLE(x,theta);
    [miuEA,sigmaEA,PEA,kEA,theta0EA] = cylinderDistOneParaMLEAppro(x,theta);
    kDiff(ng) = abs(kE-kEA);
    theta0Diff(ng) = abs(theta0E-theta0EA);
    KLDiv(ng) = cylinderDistOneParaKL(miuE,sigmaE,PE,kE,theta0E,...
        miuEA,sigmaEA,PEA,kEA,theta0EA);
    fprintf(strcat(num2str(ng),' is finished\n'));
end

figure;
plot(kV,kDiff);
set(gca,'XScale','log');
figure;
plot(kV,theta0Diff);
set(gca,'XScale','log');
figure;
plot(kV,KLDiv);
set(gca,'XScale','log');

save('D:\result-SO3Euclid\test approMLE S1 8-22-2019\varyK','miu','sigma','PV',...
    'kV','theta0','kDiff','theta0Diff','KLDiv');

%% vary sigma
sigmaV = logspace(-2,2,NGrid)';
PV = rho.*sigmaV/sqrt(k);

kDiff = zeros(NGrid,1);
theta0Diff = zeros(NGrid,1);
KLDiv = zeros(NGrid,1);
parfor ng = 1:NGrid
    [x,theta] = cylinderSampling(miu,sigmaV(ng),PV(ng,:),k,theta0,Ns);
    [miuE,sigmaE,PE,kE,theta0E] = cylinderDistOneParaMLE(x,theta);
    [miuEA,sigmaEA,PEA,kEA,theta0EA] = cylinderDistOneParaMLEAppro(x,theta);
    kDiff(ng) = abs(kE-kEA);
    theta0Diff(ng) = abs(theta0E-theta0EA);
    KLDiv(ng) = cylinderDistOneParaKL(miuE,sigmaE,PE,kE,theta0E,...
        miuEA,sigmaEA,PEA,kEA,theta0EA);
    fprintf(strcat(num2str(ng),' is finished\n'));
end

figure;
plot(sigmaV,kDiff);
set(gca,'XScale','log');
figure;
plot(sigmaV,theta0Diff);
set(gca,'XScale','log');
figure;
plot(sigmaV,KLDiv);
set(gca,'XScale','log');

save('D:\result-SO3Euclid\test approMLE S1 8-22-2019\varySigma','miu','sigmaV','PV',...
    'k','theta0','kDiff','theta0Diff','KLDiv');

%% vary P
rhoV = [zeros(NGrid,1),linspace(0,0.99,NGrid)'];
PV = rhoV*sigma/sqrt(k);

kDiff = zeros(NGrid,1);
theta0Diff = zeros(NGrid,1);
KLDiv = zeros(NGrid,1);
for ng = 1:NGrid
    [x,theta] = cylinderSampling(miu,sigma,PV(ng,:),k,theta0,Ns);
    [miuE,sigmaE,PE,kE,theta0E] = cylinderDistOneParaMLE(x,theta);
    [miuEA,sigmaEA,PEA,kEA,theta0EA] = cylinderDistOneParaMLEAppro(x,theta);
    kDiff(ng) = abs(kE-kEA);
    theta0Diff(ng) = abs(theta0E-theta0EA);
    KLDiv(ng) = cylinderDistOneParaKL(miuE,sigmaE,PE,kE,theta0E,...
        miuEA,sigmaEA,PEA,kEA,theta0EA);
    fprintf(strcat(num2str(ng),' is finished\n'));
end

figure;
plot(rhoV(:,2),kDiff);
figure;
plot(rhoV(:,2),theta0Diff);
figure;
plot(rhoV(:,2),KLDiv);

save('D:\result-SO3Euclid\test approMLE S1 8-22-2019\varyP','miu','sigma','PV',...
    'k','theta0','kDiff','theta0Diff','KLDiv');

end

