close all;
clear;

addpath('../../rotation3d');
addpath('../Generate-Path');
addpath('../Matrix-Fisher-Distribution');
addpath('../Filter');

% time parameters
T = 60;
sf = 200;

parameters.dt = 1/sf;

% noise parameters
parameters.randomWalk = 10*pi/180;
parameters.biasInstability = 500/3600*pi/180;
parameters.acceRandomWalk = 0.1;
parameters.acceBiasInstability = 200/3600;
parameters.rotMeaNoise = 0.1^2*eye(3);
parameters.posMeaNoise = 0.1^2*eye(3);
parameters.GaussMea = true;
parameters.bool_prog = false;

% generate path
[gyro,acce,RMea,pMea,RTrue,xTrue] = genTrigIMU(T,sf,parameters);

% initial parameters
parameters.RInit = RTrue(:,:,1)*expRot([0;0;pi]);
parameters.xInit = xTrue(:,1);
parameters.initRNoise = 1e10*eye(3);
parameters.initXNoise = 0.1^2*eye(12);

% filter
[RMEKF,xMEKF,G] = MEKFIMU(gyro,acce,[],pMea,parameters);
[RMFG,xMFG,MFG] = MFGAnalyticIMU(gyro,acce,[],pMea,parameters);

parameters.bool_prog = true;
[RMFGp,xMFGp,MFGp] = MFGAnalyticIMU(gyro,acce,[],pMea,parameters);

% plot error
figure; hold on;
plot(sqrt(sum(logRot(mulRot(invRot(RTrue),RMEKF),'v')).^2));
plot(sqrt(sum(logRot(mulRot(invRot(RTrue),RMFG),'v')).^2));
plot(sqrt(sum(logRot(mulRot(invRot(RTrue),RMFGp),'v')).^2));

figure; hold on;
plot(sqrt(sum((xTrue(1:3,:)-xMEKF(1:3,:)).^2)));
plot(sqrt(sum((xTrue(1:3,:)+xMFG(1:3,:)).^2)));
plot(sqrt(sum((xTrue(1:3,:)+xMFGp(1:3,:)).^2)));

figure; hold on;
plot(sqrt(sum((xTrue(4:6,:)-xMEKF(4:6,:)).^2)));
plot(sqrt(sum((xTrue(4:6,:)-xMFG(4:6,:)).^2)));
plot(sqrt(sum((xTrue(4:6,:)-xMFGp(4:6,:)).^2)));

figure; hold on;
plot(sqrt(sum((xTrue(7:9,:)-xMEKF(7:9,:)).^2)));
plot(sqrt(sum((xTrue(7:9,:)-xMFG(7:9,:)).^2)));
plot(sqrt(sum((xTrue(7:9,:)-xMFGp(7:9,:)).^2)));

figure; hold on;
plot(sqrt(sum((xTrue(10:12,:)-xMEKF(10:12,:)).^2)));
plot(sqrt(sum((xTrue(10:12,:)+xMFG(10:12,:)).^2)));
plot(sqrt(sum((xTrue(10:12,:)+xMFGp(10:12,:)).^2)));

rmpath('../../rotation3d');
rmpath('../Generate-Path');
rmpath('../Matrix-Fisher-Distribution');
rmpath('../Filter');

