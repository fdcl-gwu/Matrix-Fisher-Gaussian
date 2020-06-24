function [] = test_attMea()

addpath('..\Generate-Path');
addpath('..');
addpath('..\Matrix-Fisher-Distribution');
addpath('..\..\rotation3d');

t = 600;
sf = 150;
N = 100;

% random seeds
parfor n = 1:20
    rng(n);
end
    
parfor n = 1:N
    parameters = [];
    parameters.setting.omegaLocal = true;
    parameters.setting.GaussMea = false;
    parameters.setting.meaIsVec = false;
    parameters.setting.attMeaLocal = true;
    parameters.meaNoise = diag([200,200,200]);
    parameters.initValue.RNoise = diag([200,200,200]);
    parameters.initValue.xNoise = 0.1^2*eye(3);
    parameters.initValue.U = eye(3);
    parameters.initValue.V = eye(3);
    parameters.initValue.Miu = [0;0;0];

    [gyro,RMea,RTrue,xTrue] = genTrig(t,sf,parameters);

    [RMEKF,xMEKF,SigmaMEKF,TMEKF] = MEKF(gyro,RMea,sf,parameters);
    [RMFGA,MFGA,TMFGA] = MFGAnalytic(gyro,RMea,sf,true,parameters);
    [RMFGU,MFGU,TMFGU] = MFGUnscented(gyro,RMea,sf,true,parameters);

    parsave(n,gyro,RMea,RTrue,xTrue,parameters,RMEKF,xMEKF,SigmaMEKF,TMEKF,...
        RMFGA,MFGA,TMFGA,RMFGU,MFGU,TMFGU);

end

rmpath('..\Generate-Path');
rmpath('..');
rmpath('..\Matrix-Fisher-Distribution');
rmpath('..\..\rotation3d');

end


function [] = parsave(n,gyro,RMea,RTrue,xTrue,parameters,RMEKF,xMEKF,SigmaMEKF,TMEKF,...
    RMFGA,MFGA,TMFGA,RMFGU,MFGU,TMFGU)

save(strcat('D:\result-SO3Euclid\6-23-2020\',num2str(n)),'gyro','RMea','RTrue','xTrue','parameters','RMEKF','xMEKF','SigmaMEKF','TMEKF',...
        'RMFGA','MFGA','TMFGA','RMFGU','MFGU','TMFGU');

end

