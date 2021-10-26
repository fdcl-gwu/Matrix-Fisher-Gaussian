close all;
clear;

addpath('../../rotation3d');
addpath('../Generate-Path');
addpath('../Matrix-Fisher-Distribution');
addpath('../Filter');

N = 100;

% random seeds
if ~isempty(gcp('nocreate'))
    delete(gcp('nocreate'));
end

N_worker = 12;
parpool(N_worker);
parfor n = 1:12
    rng(n);
end

% time parameters
T = 60;
sf = 200;
sf_GPS = 20;

parameters_base.dt = 1/sf;
parameters_base.sf_GPS = sf_GPS;

% noise parameters
parameters_base.randomWalk = 1*pi/180;
parameters_base.biasInstability = 50/3600*pi/180;
parameters_base.acceRandomWalk = 0.01;
parameters_base.acceBiasInstability = 20/3600;
parameters_base.rotMeaNoise = 0.1^2*eye(3);
parameters_base.posMeaNoise = 0.1^2*eye(3);
parameters_base.GaussMea = true;
parameters_base.bool_prog = false;

% path
path = 'D:\result-SO3Euclid\bench_IMUGNSS\10-21-2021';
if ~exist(path,'dir')
    mkdir(path);
end

% pre-allocate memory
RErrorMean = zeros(N,7);
bgErrorMean = zeros(N,7);
pErrorMean = zeros(N,7);
vErrorMean = zeros(N,7);
baErrorMean = zeros(N,7);

% run
parfor i = 1:N
    % generate trajectory
    [gyro,acce,RMea,pMea,RTrue,xTrue] = genTrigIMU(T,sf,parameters_base);
    
    % initial parameters
    parameters = parameters_base;
    
    g = [0;0;9.8];
    r = cross(acce(:,1),g);
    theta = asin(sqrt(sum(r.^2))/sqrt(sum(acce(:,1).^2))/sqrt(sum(g.^2)));
    r = r/sqrt(sum(r.^2));
    parameters.RInit = expRot([0;0;pi])*expRot(theta*hat(r));

    Ns = 100000;
    aInit_rand = RTrue(:,:,1)'*g + randn(3,Ns)*parameters.acceRandomWalk*sqrt(sf);
    r = cross(aInit_rand,repmat(g,1,Ns));
    theta = asin(sqrt(sum(r.^2))./sqrt(sum(aInit_rand.^2))/sqrt(sum(g.^2)));
    r = r./sqrt(sum(r.^2));
    RInit_rand = expRot(theta.*r);
    RInit_rand = mulRot(expRot([zeros(2,Ns);(rand(1,Ns)-0.5)*2*pi]),RInit_rand);
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
    [RMEKF1,xMEKF1,G1,tMEKF1] = MEKFIMU(gyro,acce,[],pMea,parameters);
    
    parameters.initRNoise(3,3) = parameters.initRNoise(1,1);
    [RMEKF2,xMEKF2,G2,tMEKF2] = MEKFIMU(gyro,acce,[],pMea,parameters);
    
    parameters.initRNoise(3,3) = 10;
    [RMEKF3,xMEKF3,G3,tMEKF3] = MEKFIMU(gyro,acce,[],pMea,parameters);
    
    parameters.initRNoise(3,3) = 100;
    [RMEKF4,xMEKF4,G4,tMEKF4] = MEKFIMU(gyro,acce,[],pMea,parameters);
    
    parameters.initRNoise(3,3) = 1000;
    [RMEKF5,xMEKF5,G5,tMEKF5] = MEKFIMU(gyro,acce,[],pMea,parameters);
    
    parameters.GaussMea = false;
    parameters.initRNoise = FInit;
    try
        [RMFG,xMFG,MFG,tMFG] = MFGAnalyticIMU(gyro,acce,[],pMea,parameters);
    catch
        RMFG = repmat(eye(3),1,1,T*sf+1);
        xMFG = repmat(zeros(12,1),1,T*sf+1);
        MFG = [];
        tMFG = [];
    end

    parameters.bool_prog = true;
    try
        [RMFGp,xMFGp,MFGp,tMFGp] = MFGAnalyticIMU(gyro,acce,[],pMea,parameters);
    catch
        RMFGp = repmat(eye(3),1,1,T*sf+1);
        xMFGp = repmat(zeros(12,1),1,T*sf+1);
        MFGp = [];
        tMFGp = [];
    end
    
    % calculate errors
    RError = cat(1,sqrt(sum(logRot(mulRot(invRot(RTrue),RMEKF1),'v').^2)),...
        sqrt(sum(logRot(mulRot(invRot(RTrue),RMEKF2),'v').^2)),...
        sqrt(sum(logRot(mulRot(invRot(RTrue),RMEKF3),'v').^2)),...
        sqrt(sum(logRot(mulRot(invRot(RTrue),RMEKF4),'v').^2)),...
        sqrt(sum(logRot(mulRot(invRot(RTrue),RMEKF5),'v').^2)),...
        sqrt(sum(logRot(mulRot(invRot(RTrue),RMFG),'v').^2)),...
        sqrt(sum(logRot(mulRot(invRot(RTrue),RMFGp),'v').^2)));
    
    bgError = cat(1,sqrt(sum((xTrue(1:3,:)-xMEKF1(1:3,:)).^2)),...
        sqrt(sum((xTrue(1:3,:)-xMEKF2(1:3,:)).^2)),...
        sqrt(sum((xTrue(1:3,:)-xMEKF3(1:3,:)).^2)),...
        sqrt(sum((xTrue(1:3,:)-xMEKF4(1:3,:)).^2)),...
        sqrt(sum((xTrue(1:3,:)-xMEKF5(1:3,:)).^2)),...
        sqrt(sum((xTrue(1:3,:)+xMFG(1:3,:)).^2)),...
        sqrt(sum((xTrue(1:3,:)+xMFGp(1:3,:)).^2)));
    
    pError = cat(1,sqrt(sum((xTrue(4:6,:)-xMEKF1(4:6,:)).^2)),...
        sqrt(sum((xTrue(4:6,:)-xMEKF2(4:6,:)).^2)),...
        sqrt(sum((xTrue(4:6,:)-xMEKF3(4:6,:)).^2)),...
        sqrt(sum((xTrue(4:6,:)-xMEKF4(4:6,:)).^2)),...
        sqrt(sum((xTrue(4:6,:)-xMEKF5(4:6,:)).^2)),...
        sqrt(sum((xTrue(4:6,:)-xMFG(4:6,:)).^2)),...
        sqrt(sum((xTrue(4:6,:)-xMFGp(4:6,:)).^2)));
    
    vError = cat(1,sqrt(sum((xTrue(7:9,:)-xMEKF1(7:9,:)).^2)),...
        sqrt(sum((xTrue(7:9,:)-xMEKF2(7:9,:)).^2)),...
        sqrt(sum((xTrue(7:9,:)-xMEKF3(7:9,:)).^2)),...
        sqrt(sum((xTrue(7:9,:)-xMEKF4(7:9,:)).^2)),...
        sqrt(sum((xTrue(7:9,:)-xMEKF5(7:9,:)).^2)),...
        sqrt(sum((xTrue(7:9,:)-xMFG(7:9,:)).^2)),...
        sqrt(sum((xTrue(7:9,:)-xMFGp(7:9,:)).^2)));
    
    baError = cat(1,sqrt(sum((xTrue(10:12,:)-xMEKF1(10:12,:)).^2)),...
        sqrt(sum((xTrue(10:12,:)-xMEKF2(10:12,:)).^2)),...
        sqrt(sum((xTrue(10:12,:)-xMEKF3(10:12,:)).^2)),...
        sqrt(sum((xTrue(10:12,:)-xMEKF4(10:12,:)).^2)),...
        sqrt(sum((xTrue(10:12,:)-xMEKF5(10:12,:)).^2)),...
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
    parsave(filepath,gyro,acce,RMea,pMea,RTrue,xTrue,RMEKF1,xMEKF1,G1,tMEKF1,...
        RMEKF2,xMEKF2,G2,tMEKF2,RMEKF3,xMEKF3,G3,tMEKF3,...
        RMEKF4,xMEKF4,G4,tMEKF4,RMEKF5,xMEKF5,G5,tMEKF5,...
        RMFG,xMFG,MFG,tMFG,RMFGp,xMFGp,MFGp,tMFGp,...
        RError,bgError,pError,vError,baError);
end

% save error
filepath = strcat(path,'\error.mat');
parsave(filepath,RErrorMean,bgErrorMean,pErrorMean,vErrorMean,baErrorMean);

% barplot for errors
figure; hold on;
bar(mean(RErrorMean));
errorbar(mean(RErrorMean),std(RErrorMean),'k','LineStyle','none');
title(strcat('Attitude Error (N=',num2str(N),')'));
xticks(1:7);
xticklabels({'MEKF1','MEKF2','MEKF3','MEKF4','MEKF5','MFG','MFGp'});
ylabel('radian');

figure; hold on;
bar(mean(bgErrorMean));
errorbar(mean(bgErrorMean),std(bgErrorMean),'k','LineStyle','none');
title(strcat('Gyroscope Bias Error (N=',num2str(N),')'));
xticks(1:7);
xticklabels({'MEKF1','MEKF2','MEKF3','MEKF4','MEKF5','MFG','MFGp'});
ylabel('radian/s');

figure; hold on;
bar(mean(pErrorMean));
errorbar(mean(pErrorMean),std(pErrorMean),'k','LineStyle','none');
title(strcat('Position Error (N=',num2str(N),')'));
xticks(1:7);
xticklabels({'MEKF1','MEKF2','MEKF3','MEKF4','MEKF5','MFG','MFGp'});
ylabel('meter');

figure; hold on;
bar(mean(vErrorMean));
errorbar(mean(vErrorMean),std(vErrorMean),'k','LineStyle','none');
title(strcat('Velocity Error (N=',num2str(N),')'));
xticks(1:7);
xticklabels({'MEKF1','MEKF2','MEKF3','MEKF4','MEKF5','MFG','MFGp'});
ylabel('meter/s');

figure; hold on;
bar(mean(baErrorMean));
errorbar(mean(baErrorMean),std(baErrorMean),'k','LineStyle','none');
title(strcat('Accelerometer Bias Error (N=',num2str(N),')'));
xticks(1:7);
xticklabels({'MEKF1','MEKF2','MEKF3','MEKF4','MEKF5','MFG','MFGp'});
ylabel('meter/s^2');

rmpath('../../rotation3d');
rmpath('../Generate-Path');
rmpath('../Matrix-Fisher-Distribution');
rmpath('../Filter');
