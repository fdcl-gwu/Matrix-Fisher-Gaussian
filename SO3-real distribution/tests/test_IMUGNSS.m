close all;
clear;

addpath('../../rotation3d');
addpath('../Generate-Path');
addpath('../Matrix-Fisher-Distribution');
addpath('../Filter');

% rng(0);

% time parameters
T = 300;
sf = 100;
sf_GPS = 10;

parameters.dt = 1/sf;
parameters.sf_GPS = sf_GPS;

% noise parameters
parameters.randomWalk = 1*pi/180;
parameters.biasInstability = 50/3600*pi/180;
parameters.acceRandomWalk = 0.01;
parameters.acceBiasInstability = 20/3600;
parameters.acceDynamics = 0.1^2*sqrt(sf);
parameters.rotMeaNoise = 0.1^2*eye(3);
parameters.posMeaNoise = 1^2*eye(3);
parameters.gravMeaNoise = 5^2*eye(3);
parameters.GaussMea = true;
parameters.noGPS = [200,15000];
parameters.useGrav = false;
parameters.bool_prog = true;

% generate path
[gyro,acce,RMea,pMea,RTrue,xTrue] = genTrigIMU(T,sf,parameters);

% initial parameters
% g = [0;0;9.8];
% r = cross(acce(:,1),g);
% theta = asin(sqrt(sum(r.^2))/sqrt(sum(acce(:,1).^2))/sqrt(sum(g.^2)));
% r = r/sqrt(sum(r.^2));
% parameters.RInit = expRot([0;0;pi])*expRot(theta*hat(r));
% 
% N = 100000;
% aInit_rand = RTrue(:,:,1)'*g + randn(3,N)*parameters.acceRandomWalk*sqrt(sf);
% r = cross(aInit_rand,repmat(g,1,N));
% theta = asin(sqrt(sum(r.^2))./sqrt(sum(aInit_rand.^2))/sqrt(sum(g.^2)));
% r = r./sqrt(sum(r.^2));
% RInit_rand = expRot(theta.*r);
% RInit_rand = mulRot(expRot([zeros(2,N);(rand(1,N)-0.5)*2*pi]),RInit_rand);
% RInit_error = mulRot(parameters.RInit',RInit_rand);
% 
% ERInit_error = mean(RInit_error,3);
% [UInit,DInit,VInit] = psvd(ERInit_error);
% SInit = diag(pdf_MF_M2S(diag(DInit)));
% FInit = parameters.RInit*(UInit*VInit')'*UInit*SInit*VInit';
% 
% rv = logRot(RInit_error,'v');
% SigmaRInit = cov(rv');
% sigma = eig(SigmaRInit);
% parameters.initRNoise = diag([mean(sigma(1:2)),mean(sigma(1:2)),sigma(3)]);
% 
% parameters.xInit = [zeros(3,1);pMea(:,1);zeros(6,1)];
% parameters.initXNoise = [0.01^2*eye(3),zeros(3,9);
%     zeros(3,3),parameters.posMeaNoise,zeros(3,6);zeros(6,6),0.01^2*eye(6)];

parameters.RInit = RTrue(:,:,1)*expRot(0.1*randn(3,1));
parameters.initRNoise = 0.1^2*eye(3);

parameters.xInit = [zeros(3,1);pMea(:,1);zeros(6,1)];
parameters.initXNoise = [0.01^2*eye(3),zeros(3,9);
    zeros(3,3),parameters.posMeaNoise,zeros(3,6);zeros(6,6),0.01^2*eye(6)];

% filter
[RMEKF,xMEKF,G] = MEKFIMU(gyro,acce,[],pMea,parameters);

parameters.GaussMea = false;
parameters.initRNoise = 50*eye(3)*parameters.RInit;
[RMFGI,xMFGI,MFGI] = MFGAnalyticIMU(gyro,acce,[],pMea,true,parameters);
[RMFGB,xMFGB,MFGB] = MFGAnalyticIMU(gyro,acce,[],pMea,false,parameters);

% plot error
figure; hold on;
plot(sqrt(sum(logRot(mulRot(invRot(RTrue),RMEKF),'v').^2)));
plot(sqrt(sum(logRot(mulRot(invRot(RTrue),RMFGI),'v').^2)));
plot(sqrt(sum(logRot(mulRot(invRot(RTrue),RMFGB),'v').^2)));

figure; hold on;
plot(sqrt(sum((xTrue(1:3,:)-xMEKF(1:3,:)).^2)));
plot(sqrt(sum((xTrue(1:3,:)-xMFGI(1:3,:)).^2)));
plot(sqrt(sum((xTrue(1:3,:)-xMFGB(1:3,:)).^2)));

figure; hold on;
plot(sqrt(sum((xTrue(4:6,:)-xMEKF(4:6,:)).^2)));
plot(sqrt(sum((xTrue(4:6,:)-xMFGI(4:6,:)).^2)));
plot(sqrt(sum((xTrue(4:6,:)-xMFGB(4:6,:)).^2)));

figure; hold on;
plot(sqrt(sum((xTrue(7:9,:)-xMEKF(7:9,:)).^2)));
plot(sqrt(sum((xTrue(7:9,:)-xMFGI(7:9,:)).^2)));
plot(sqrt(sum((xTrue(7:9,:)-xMFGB(7:9,:)).^2)));

figure; hold on;
plot(sqrt(sum((xTrue(10:12,:)-xMEKF(10:12,:)).^2)));
plot(sqrt(sum((xTrue(10:12,:)-xMFGI(10:12,:)).^2)));
plot(sqrt(sum((xTrue(10:12,:)-xMFGB(10:12,:)).^2)));

rmpath('../../rotation3d');
rmpath('../Generate-Path');
rmpath('../Matrix-Fisher-Distribution');
rmpath('../Filter');

