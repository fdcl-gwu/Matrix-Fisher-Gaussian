function [] = test()

addpath('..\Generate-Path');

rng(0.01);
N = 60;

t = 60;
sf = 50;

parfor n = 1:N
    [gyro,RMea,RTrue,xTrue] = genTrigWithBias(t,sf);
    [RMEKF,xMEKF,G] = MEKFWithBias(gyro,RTrue(:,:,1),RMea,sf);
    [RMFG,MFG] = MFGUnscentedWithBias(gyro,RTrue(:,:,1),RMea,sf);
    
    parsave(n,gyro,RMea,RTrue,xTrue,RMEKF,xMEKF,G,RMFG,MFG);
end

rmpath('..\Generate-Path');

end


function [] = parsave(n,gyro,RMea,RTrue,xTrue,RMEKF,xMEKF,G,RMFG,MFG)

save(strcat('D:\result-SO3Euclid\Unscented 9-10-2019\',num2str(n),'.mat'),...
    'gyro','RMea','RTrue','xTrue','RMEKF','xMEKF','G','RMFG','MFG','-v7.3');

end

