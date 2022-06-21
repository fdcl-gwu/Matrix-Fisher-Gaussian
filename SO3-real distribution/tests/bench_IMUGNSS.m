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

N_worker = 16;
parpool(N_worker);
parfor n = 1:16
    rng(n);
end

% time parameters
T = 60;
sf = 200;
sf_GPS = 10;

parameters_base.dt = 1/sf;
parameters_base.sf_GPS = sf_GPS;

% noise parameters
parameters_base.randomWalk = 1*pi/180;
parameters_base.biasInstability = 50/3600*pi/180;
parameters_base.acceRandomWalk = 0.01;
parameters_base.acceBiasInstability = 20/3600;
parameters_base.rotMeaNoise = 0.1^2*eye(3);
parameters_base.posMeaNoise = 1^2*eye(3);
parameters_base.gravMeaNoise = 6^2*eye(3);
parameters_base.GaussMea = true;
parameters_base.useGrav = false;
parameters_base.bool_prog = true;

% path
path = 'D:\result-SO3Euclid\bench_IMUGNSS\6-21-2022';
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

    parameters.xInit = [zeros(3,1);pMea(:,1);zeros(6,1)];
    parameters.initXNoise = [0.01^2*eye(3),zeros(3,9);
        zeros(3,3),parameters.posMeaNoise,zeros(3,6);zeros(6,6),0.01^2*eye(6)];
    
    % filter
    parameters.initRNoise(3,3) = 1000;
    [RMEKF,xMEKF,G,tMEKF] = MEKFIMU(gyro,acce,[],pMea,parameters);
    
    parameters.GaussMea = false;
    parameters.initRNoise = FInit;
    try
        [RMFGI,xMFGI,MFGI,tMFGI] = MFGAnalyticIMU(gyro,acce,[],pMea,true,parameters);
    catch
        RMFGI = repmat(eye(3),1,1,T*sf+1);
        xMFGI = repmat(zeros(12,1),1,T*sf+1);
        MFGI = [];
        tMFGI = [];
    end
    
    try
        [RMFGB,xMFGB,MFGB,tMFGB] = MFGAnalyticIMU(gyro,acce,[],pMea,false,parameters);
    catch
        RMFGB = repmat(eye(3),1,1,T*sf+1);
        xMFGB = repmat(zeros(12,1),1,T*sf+1);
        MFGB = [];
        tMFGB = [];
    end
    
    % calculate errors
    RError = cat(1,sqrt(sum(logRot(mulRot(invRot(RTrue),RMEKF),'v').^2)),...
        sqrt(sum(logRot(mulRot(invRot(RTrue),RMFGI),'v').^2)),...
        sqrt(sum(logRot(mulRot(invRot(RTrue),RMFGB),'v').^2)));
    
    bgError = cat(1,sqrt(sum((xTrue(1:3,:)-xMEKF(1:3,:)).^2)),...
        sqrt(sum((xTrue(1:3,:)+xMFGI(1:3,:)).^2)),...
        sqrt(sum((xTrue(1:3,:)+xMFGB(1:3,:)).^2)));
    
    pError = cat(1,sqrt(sum((xTrue(4:6,:)-xMEKF(4:6,:)).^2)),...
        sqrt(sum((xTrue(4:6,:)-xMFGI(4:6,:)).^2)),...
        sqrt(sum((xTrue(4:6,:)-xMFGB(4:6,:)).^2)));
    
    vError = cat(1,sqrt(sum((xTrue(7:9,:)-xMEKF(7:9,:)).^2)),...
        sqrt(sum((xTrue(7:9,:)-xMFGI(7:9,:)).^2)),...
        sqrt(sum((xTrue(7:9,:)-xMFGB(7:9,:)).^2)));
    
    baError = cat(1,sqrt(sum((xTrue(10:12,:)-xMEKF(10:12,:)).^2)),...
        sqrt(sum((xTrue(10:12,:)+xMFGI(10:12,:)).^2)),...
        sqrt(sum((xTrue(10:12,:)+xMFGB(10:12,:)).^2)));
    
    % mean error
    RErrorMean(i,:) = mean(RError,2)';
    bgErrorMean(i,:) = mean(bgError,2)';
    pErrorMean(i,:) = mean(pError,2)';
    vErrorMean(i,:) = mean(vError,2)';
    baErrorMean(i,:) = mean(baError,2)';
    
    % save data
    filepath = strcat(path,'\',num2str(i),'.mat');
    parsave(filepath,gyro,acce,RMea,pMea,RTrue,xTrue,RMEKF,xMEKF,G,tMEKF,...
        RMFGI,xMFGI,MFGI,tMFGI,RMFGB,xMFGB,MFGB,tMFGB,...
        RError,bgError,pError,vError,baError,parameters);
end

% save error
filepath = strcat(path,'\error.mat');
parsave(filepath,RErrorMean,bgErrorMean,pErrorMean,vErrorMean,baErrorMean);

rmpath('../../rotation3d');
rmpath('../Generate-Path');
rmpath('../Matrix-Fisher-Distribution');
rmpath('../Filter');
