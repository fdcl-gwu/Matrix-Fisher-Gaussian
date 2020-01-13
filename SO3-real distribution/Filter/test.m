function [] = test()

addpath('..\Generate-Path');

rng(0.01);
N = 60;

t = 60;
sf = 50;

parfor n = 1:N
    % general settings
    parameters = [];
    parameters.t = 60;
    parameters.dt = 1/sf;
    parameters.randomWalk = 10*pi/180;
    parameters.biasInstability = 500/3600*pi/180;
    parameters.rotMeaNoise = 0.2;
    parameters.initRsigma = 0.2;
    parameters.initXsigma = 0.01;
    parameters.RInit = eye(3);
    parameters.xInit = [0;0;0];
    
    [gyro,RMea,RTrue,xTrue] = genTrigWithBias(t,sf,parameters);
    
    % stachastic settings
    parameters.RInit = RMea(:,:,1);
    parameters.xInit = randn(3,1)*parameters.initXsigma;
    
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

save(strcat('D:\result-SO3Euclid\1-11-2020\',num2str(n),'.mat'),...
    'gyro','RMea','RTrue','xTrue','parameters','RMEKF','xMEKF','G','TMEKF',...
    'RMFGA','MFGA','TMFGA','RMFGU','MFGU','TMFGU','-v7.3');

end

