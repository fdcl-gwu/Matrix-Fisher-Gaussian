function [] = test()

addpath('..\Generate-Path');

parpool(15);
parfor n = 1:15
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
    parameters.rotMeaNoise = 20;
    parameters.initRNoise = 20;
    parameters.initXNoise = 0.1;
    parameters.RInit = eye(3);
    parameters.xInit = [0;0;0];
    
    [gyro,RMea,RTrue,xTrue] = genTrigWithBias(t,sf,parameters);
    
    % stachastic settings
    parameters.RInit = RTrue(:,:,1)*expRot([pi,0,0]);
    parameters.xInit = [0.2;0.2;0.2];
    
    [RMEKF,xMEKF,G,TMEKF] = MEKFWithBias(gyro,RMea,parameters);
    [RMFGA,MFGA,TMFGA] = MFGAnaWithBias(gyro,RMea,parameters);
    [RMFGU,MFGU,TMFGU] = MFGUnscentedWithBias(gyro,RMea,parameters);
    
    parsave(n,gyro,RMea,RTrue,xTrue,parameters,RMEKF,xMEKF,G,TMEKF,...
        RMFGA,MFGA,TMFGA,RMFGU,MFGU,TMFGU);
end

rmpath('..\Generate-Path');

end


function [] = parsave(n,gyro,RMea,RTrue,xTrue,parameters,RMEKF,xMEKF,G,TMEKF,...
    RMFGA,MFGA,TMFGA,RMFGU,MFGU,TMFGU)

save(strcat('D:\result-SO3Euclid\1-16-2020-3\',num2str(n),'.mat'),...
    'gyro','RMea','RTrue','xTrue','parameters','RMEKF','xMEKF','G','TMEKF',...
    'RMFGA','MFGA','TMFGA','RMFGU','MFGU','TMFGU','-v7.3');

end

