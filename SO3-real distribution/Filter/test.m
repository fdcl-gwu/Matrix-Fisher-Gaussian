function [] = test()

addpath('..\Generate-Path');
addpath('..');
addpath('..\Matrix-Fisher-Distribution');
addpath('..\..\rotation3d');

t = 60;
sf = 150;
N = 100;

% random seeds
parfor n = 1:20
    rng(n);
end
    
parfor n = 1:N
    parameters = [];
    parameters.setting.omegaLocal = true;
    parameters.setting.gyroFail = true;
    parameters.setting.GaussMea = false;
    parameters.setting.meaIsVec = true;
    parameters.setting.vecRefInertial = false;
    parameters.setting.nVecRef = 1;
    parameters.setting.vRef = [1;0;0];
    parameters.meaNoise = 200;
    parameters.initValue.RNoise = diag([200,200,200]);
    parameters.initValue.xNoise = 0.1^2*eye(3);
    parameters.initValue.U = eye(3);
    parameters.initValue.V = eye(3);
    parameters.initValue.Miu = [0;0;0];

    [gyro,RMea,RTrue,xTrue] = genTrig(t,sf,parameters);

    [RMEKF,xMEKF,SigmaMEKF,TMEKF] = MEKF(gyro,RMea,sf,parameters);
    [RMFGAQS,MFGAQS,TMFGAQS] = MFGAnalytic(gyro,RMea,sf,true,parameters);
    [RMFGUQS,MFGUQS,TMFGUQS] = MFGUnscented(gyro,RMea,sf,true,parameters);
    [RMFGASQ,MFGASQ,TMFGASQ] = MFGAnalytic(gyro,RMea,sf,false,parameters);
    [RMFGUSQ,MFGUSQ,TMFGUSQ] = MFGUnscented(gyro,RMea,sf,false,parameters);

    parsave(n,gyro,RMea,RTrue,xTrue,parameters,RMEKF,xMEKF,SigmaMEKF,TMEKF,...
        RMFGAQS,MFGAQS,TMFGAQS,RMFGUQS,MFGUQS,TMFGUQS,...
        RMFGASQ,MFGASQ,TMFGASQ,RMFGUSQ,MFGUSQ,TMFGUSQ);

end

rmpath('..\Generate-Path');
rmpath('..');
rmpath('..\Matrix-Fisher-Distribution');
rmpath('..\..\rotation3d');

end


function [] = parsave(n,gyro,RMea,RTrue,xTrue,parameters,RMEKF,xMEKF,SigmaMEKF,TMEKF,...
    RMFGAQS,MFGAQS,TMFGAQS,RMFGUQS,MFGUQS,TMFGUQS,...
    RMFGASQ,MFGASQ,TMFGASQ,RMFGUSQ,MFGUSQ,TMFGUSQ)

save(strcat('D:\result-SO3Euclid\7-14-2020-4\',num2str(n)),'gyro','RMea','RTrue','xTrue','parameters','RMEKF','xMEKF','SigmaMEKF','TMEKF',...
        'RMFGAQS','MFGAQS','TMFGAQS','RMFGUQS','MFGUQS','TMFGUQS',...
        'RMFGASQ','MFGASQ','TMFGASQ','RMFGUSQ','MFGUSQ','TMFGUSQ');

end

