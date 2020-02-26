function [] = test()

addpath('..\Generate-Path');

NWorker = 20;
parpool(NWorker);
parfor n = 1:NWorker
    rng(n);
end

N = 60;

t = 60;
sf = 150;

parfor n = 1:N
    % general settings
    parameters = [];
    parameters.t = 60;
    parameters.dt = 1/sf;
    parameters.randomWalk = 10*pi/180;
    parameters.biasInstability = 500/3600*pi/180;
    parameters.GaussMea = false;
    parameters.rotMeaNoise = 12*eye(3);
    parameters.initRNoise = 0*eye(3);
    parameters.initXNoise = 0.2^2*eye(3);
    parameters.RInit = eye(3);
    parameters.xInit = [0;0;0];
    
    [gyro,RMea,RTrue,xTrue] = genTrigWithBias(t,sf,parameters);
    
    % stachastic settings
    parameters.RInit = RTrue(:,:,1)*expRot([pi,0,0]);
    parameters.xInit = [0.2;0.2;0.2];
    
    [RMEKF,xMEKF,G,TMEKF] = MEKF(gyro,RMea,parameters);
    [RMFGA,MFGA,TMFGA] = MFGAnalytic(gyro,RMea,parameters);
    [RMFGU,MFGU,TMFGU] = MFGUnscented(gyro,RMea,parameters);
    
    parsave(n,gyro,RMea,RTrue,xTrue,parameters,RMEKF,xMEKF,G,TMEKF,...
        RMFGA,MFGA,TMFGA,RMFGU,MFGU,TMFGU);
end

rmpath('..\Generate-Path');

end


function [] = parsave(n,gyro,RMea,RTrue,xTrue,parameters,RMEKF,xMEKF,G,TMEKF,...
    RMFGA,MFGA,TMFGA,RMFGU,MFGU,TMFGU)

save(strcat('D:\result-SO3Euclid\2-26-2020\',num2str(n),'.mat'),...
    'gyro','RMea','RTrue','xTrue','parameters','RMEKF','xMEKF','G','TMEKF',...
    'RMFGA','MFGA','TMFGA','RMFGU','MFGU','TMFGU','-v7.3');

end

