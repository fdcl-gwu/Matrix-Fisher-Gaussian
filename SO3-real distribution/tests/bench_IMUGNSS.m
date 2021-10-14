close all;
clear;

addpath('../../rotation3d');
addpath('../Generate-Path');
addpath('../Matrix-Fisher-Distribution');
addpath('../Filter');

rng(1);
N = 2;

% time parameters
T = 60;
sf = 200;

parameters_base.dt = 1/sf;

% noise parameters
parameters_base.randomWalk = 10*pi/180;
parameters_base.biasInstability = 500/3600*pi/180;
parameters_base.acceRandomWalk = 0.1;
parameters_base.acceBiasInstability = 200/3600;
parameters_base.rotMeaNoise = 0.1^2*eye(3);
parameters_base.posMeaNoise = 0.1^2*eye(3);
parameters_base.GaussMea = true;
parameters_base.bool_prog = false;

% path
path = 'D:\result-SO3Euclid\bench_IMUGNSS\10-14-2021';
if ~exist(path,'dir')
    mkdir(path);
end

% pre-allocate memory
RErrorMean = zeros(N,3);
bgErrorMean = zeros(N,3);
pErrorMean = zeros(N,3);
vErrorMean = zeros(N,3);
baErrorMean = zeros(N,3);

% run
parfor i = 1:N
    % generate trajectory
    [gyro,acce,RMea,pMea,RTrue,xTrue] = genTrigIMU(T,sf,parameters_base);
    
    % initial parameters
    parameters = parameters_base;
    parameters.RInit = RTrue(:,:,1)*expRot([0;0;pi]);
    parameters.xInit = xTrue(:,1);
    parameters.initRNoise = 1e10*eye(3);
    parameters.initXNoise = 0.1^2*eye(12);
    
    % filter
    [RMEKF,xMEKF,G,tMEKF] = MEKFIMU(gyro,acce,[],pMea,parameters);
    [RMFG,xMFG,MFG,tMFG] = MFGAnalyticIMU(gyro,acce,[],pMea,parameters);

    parameters.bool_prog = true;
    [RMFGp,xMFGp,MFGp,tMFGp] = MFGAnalyticIMU(gyro,acce,[],pMea,parameters);
    
    % calculate errors
    RError = cat(1,sqrt(sum(logRot(mulRot(invRot(RTrue),RMEKF),'v').^2)),...
        sqrt(sum(logRot(mulRot(invRot(RTrue),RMFG),'v').^2)),...
        sqrt(sum(logRot(mulRot(invRot(RTrue),RMFGp),'v').^2)));
    
    bgError = cat(1,sqrt(sum((xTrue(1:3,:)-xMEKF(1:3,:)).^2)),...
        sqrt(sum((xTrue(1:3,:)+xMFG(1:3,:)).^2)),...
        sqrt(sum((xTrue(1:3,:)+xMFGp(1:3,:)).^2)));
    
    pError = cat(1,sqrt(sum((xTrue(4:6,:)-xMEKF(4:6,:)).^2)),...
        sqrt(sum((xTrue(4:6,:)-xMFG(4:6,:)).^2)),...
        sqrt(sum((xTrue(4:6,:)-xMFGp(4:6,:)).^2)));
    
    vError = cat(1,sqrt(sum((xTrue(7:9,:)-xMEKF(7:9,:)).^2)),...
        sqrt(sum((xTrue(7:9,:)-xMFG(7:9,:)).^2)),...
        sqrt(sum((xTrue(7:9,:)-xMFGp(7:9,:)).^2)));
    
    baError = cat(1,sqrt(sum((xTrue(10:12,:)-xMEKF(10:12,:)).^2)),...
        sqrt(sum((xTrue(10:12,:)+xMFG(10:12,:)).^2)),...
        sqrt(sum((xTrue(10:12,:)+xMFGp(10:12,:)).^2)));
    
    % mean error
    RErrorMean(i,:) = mean(RError,2)';
    bgErrorMean(i,:) = mean(bgError,2)';
    pErrorMean(i,:) = mean(pError,2)';
    vErrorMean(i,:) = mean(vError,2)';
    baErrorMean(i,:) = mean(baError,2)';
    
    % figure;
    figure;
    plot(RError');
    
    figure;
    plot(bgError');
    
    figure;
    plot(pError');
    
    figure;
    plot(vError');
    
    figure;
    plot(baError');
    
    % save data
    filepath = strcat(path,'\',num2str(i),'.mat');
    parsave(filepath,gyro,acce,RMea,pMea,RTrue,xTrue,RMEKF,xMEKF,G,tMEKF,...
        RMFG,xMFG,MFG,tMFG,RMFGp,xMFGp,MFGp,tMFGp,...
        RError,bgError,pError,vError,baError);
end

% barplot for errors
figure; hold on;
bar(mean(RErrorMean));
errorbar(mean(RErrorMean),std(RErrorMean),'k','LineStyle','none');
title(strcat('Attitude Error (N=',num2str(N),')'));
xticks([1,2,3]);
xticklabels({'Gauss','MFG','MFGp'});
ylabel('radian');

figure; hold on;
bar(mean(bgErrorMean));
errorbar(mean(bgErrorMean),std(bgErrorMean),'k','LineStyle','none');
title(strcat('Gyroscope Bias Error (N=',num2str(N),')'));
xticks([1,2,3]);
xticklabels({'Gauss','MFG','MFGp'});
ylabel('radian/s');

figure; hold on;
bar(mean(pErrorMean));
errorbar(mean(pErrorMean),std(pErrorMean),'k','LineStyle','none');
title(strcat('Position Error (N=',num2str(N),')'));
xticks([1,2,3]);
xticklabels({'Gauss','MFG','MFGp'});
ylabel('meter');

figure; hold on;
bar(mean(vErrorMean));
errorbar(mean(vErrorMean),std(vErrorMean),'k','LineStyle','none');
title(strcat('Velocity Error (N=',num2str(N),')'));
xticks([1,2,3]);
xticklabels({'Gauss','MFG','MFGp'});
ylabel('meter/s');

figure; hold on;
bar(mean(baErrorMean));
errorbar(mean(baErrorMean),std(baErrorMean),'k','LineStyle','none');
title(strcat('Accelerometer Bias Error (N=',num2str(N),')'));
xticks([1,2,3]);
xticklabels({'Gauss','MFG','MFGp'});
ylabel('meter/s^2');

rmpath('../../rotation3d');
rmpath('../Generate-Path');
rmpath('../Matrix-Fisher-Distribution');
rmpath('../Filter');
