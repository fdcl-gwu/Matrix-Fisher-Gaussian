close all;
clear;

addpath('../../rotation3d');
addpath('../Generate-Path');
addpath('../Matrix-Fisher-Distribution');
addpath('../Filter');

% time parameters
T = 60;
sf = 200;
sf_GPS = 40;

parameters.dt = 1/sf;
parameters.sf_GPS = sf_GPS;

% noise parameters
parameters.randomWalk = 1*pi/180;
parameters.biasInstability = 50/3600*pi/180;
parameters.acceRandomWalk = 0.01;
parameters.acceBiasInstability = 20/3600;
parameters.rotMeaNoise = 0.1^2*eye(3);
parameters.posMeaNoise = 0.1^2*eye(3);
parameters.GaussMea = true;
parameters.bool_prog = false;

% generate path
[gyro,acce,RMea,pMea,RTrue,xTrue] = genTrigIMU(T,sf,parameters);

% initial parameters
g = [0;0;9.8];
r = cross(acce(:,1),g);
theta = asin(sqrt(sum(r.^2))/sqrt(sum(acce(:,1).^2))/sqrt(sum(g.^2)));
r = r/sqrt(sum(r.^2));
parameters.RInit = expRot([0;0;pi])*expRot(theta*hat(r));

N = 100000;
aInit_rand = RTrue(:,:,1)'*g + randn(3,N)*parameters.acceRandomWalk*sqrt(sf);
r = cross(aInit_rand,repmat(g,1,N));
theta = asin(sqrt(sum(r.^2))./sqrt(sum(aInit_rand.^2))/sqrt(sum(g.^2)));
r = r./sqrt(sum(r.^2));
RInit_rand = expRot(theta.*r);
RInit_rand = mulRot(expRot([zeros(2,N);(rand(1,N)-0.5)*2*pi]),RInit_rand);
RInit_error = mulRot(parameters.RInit',RInit_rand);

ERInit_error = mean(RInit_error,3);
[UInit,DInit,VInit] = psvd(ERInit_error);
SInit = diag(pdf_MF_M2S(diag(DInit)));
FInit = parameters.RInit*(UInit*VInit')'*UInit*SInit*VInit';

rv = logRot(RInit_error,'v');
SigmaRInit = cov(rv');
sigma = eig(SigmaRInit);
parameters.initRNoise = diag([mean(sigma(1:2)),mean(sigma(1:2)),sigma(3)]);

parameters.xInit = [xTrue(1:3,1);pMea(:,1);xTrue(7:12,1)];
parameters.initXNoise = [0.01^2*eye(3),zeros(3,9);
    zeros(3,3),parameters.posMeaNoise,zeros(3,6);zeros(6,6),0.01^2*eye(6)];

% filter
[RMEKF1,xMEKF1,G1] = MEKFIMU(gyro,acce,[],pMea,parameters);

parameters.initRNoise(3,3) = parameters.initRNoise(1,1);
[RMEKF2,xMEKF2,G2] = MEKFIMU(gyro,acce,[],pMea,parameters);

parameters.initRNoise(3,3) = 10;
[RMEKF3,xMEKF3,G3] = MEKFIMU(gyro,acce,[],pMea,parameters);

parameters.initRNoise(3,3) = 100;
[RMEKF4,xMEKF4,G4] = MEKFIMU(gyro,acce,[],pMea,parameters);

parameters.initRNoise(3,3) = 1000;
[RMEKF5,xMEKF5,G5] = MEKFIMU(gyro,acce,[],pMea,parameters);

parameters.GaussMea = false;
parameters.initRNoise = FInit;
[RMFG,xMFG,MFG] = MFGAnalyticIMU(gyro,acce,[],pMea,parameters);

parameters.bool_prog = true;
[RMFGp,xMFGp,MFGp] = MFGAnalyticIMU(gyro,acce,[],pMea,parameters);

% plot error
figure; hold on;
plot(sqrt(sum(logRot(mulRot(invRot(RTrue),RMEKF1),'v')).^2));
plot(sqrt(sum(logRot(mulRot(invRot(RTrue),RMEKF2),'v')).^2));
plot(sqrt(sum(logRot(mulRot(invRot(RTrue),RMEKF3),'v')).^2));
plot(sqrt(sum(logRot(mulRot(invRot(RTrue),RMEKF4),'v')).^2));
plot(sqrt(sum(logRot(mulRot(invRot(RTrue),RMEKF5),'v')).^2));
plot(sqrt(sum(logRot(mulRot(invRot(RTrue),RMFG),'v')).^2));
plot(sqrt(sum(logRot(mulRot(invRot(RTrue),RMFGp),'v')).^2));

figure; hold on;
plot(sqrt(sum((xTrue(1:3,:)-xMEKF1(1:3,:)).^2)));
plot(sqrt(sum((xTrue(1:3,:)-xMEKF2(1:3,:)).^2)));
plot(sqrt(sum((xTrue(1:3,:)-xMEKF3(1:3,:)).^2)));
plot(sqrt(sum((xTrue(1:3,:)-xMEKF4(1:3,:)).^2)));
plot(sqrt(sum((xTrue(1:3,:)-xMEKF5(1:3,:)).^2)));
plot(sqrt(sum((xTrue(1:3,:)+xMFG(1:3,:)).^2)));
plot(sqrt(sum((xTrue(1:3,:)+xMFGp(1:3,:)).^2)));

figure; hold on;
plot(sqrt(sum((xTrue(4:6,:)-xMEKF1(4:6,:)).^2)));
plot(sqrt(sum((xTrue(4:6,:)-xMEKF2(4:6,:)).^2)));
plot(sqrt(sum((xTrue(4:6,:)-xMEKF3(4:6,:)).^2)));
plot(sqrt(sum((xTrue(4:6,:)-xMEKF4(4:6,:)).^2)));
plot(sqrt(sum((xTrue(4:6,:)-xMEKF5(4:6,:)).^2)));
plot(sqrt(sum((xTrue(4:6,:)-xMFG(4:6,:)).^2)));
plot(sqrt(sum((xTrue(4:6,:)-xMFGp(4:6,:)).^2)));

figure; hold on;
plot(sqrt(sum((xTrue(7:9,:)-xMEKF1(7:9,:)).^2)));
plot(sqrt(sum((xTrue(7:9,:)-xMEKF2(7:9,:)).^2)));
plot(sqrt(sum((xTrue(7:9,:)-xMEKF3(7:9,:)).^2)));
plot(sqrt(sum((xTrue(7:9,:)-xMEKF4(7:9,:)).^2)));
plot(sqrt(sum((xTrue(7:9,:)-xMEKF5(7:9,:)).^2)));
plot(sqrt(sum((xTrue(7:9,:)-xMFG(7:9,:)).^2)));
plot(sqrt(sum((xTrue(7:9,:)-xMFGp(7:9,:)).^2)));

figure; hold on;
plot(sqrt(sum((xTrue(10:12,:)-xMEKF1(10:12,:)).^2)));
plot(sqrt(sum((xTrue(10:12,:)-xMEKF2(10:12,:)).^2)));
plot(sqrt(sum((xTrue(10:12,:)-xMEKF3(10:12,:)).^2)));
plot(sqrt(sum((xTrue(10:12,:)-xMEKF4(10:12,:)).^2)));
plot(sqrt(sum((xTrue(10:12,:)-xMEKF5(10:12,:)).^2)));
plot(sqrt(sum((xTrue(10:12,:)+xMFG(10:12,:)).^2)));
plot(sqrt(sum((xTrue(10:12,:)+xMFGp(10:12,:)).^2)));

rmpath('../../rotation3d');
rmpath('../Generate-Path');
rmpath('../Matrix-Fisher-Distribution');
rmpath('../Filter');

