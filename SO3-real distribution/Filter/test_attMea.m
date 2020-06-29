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
    parameters.setting.omegaLocal = false;
    parameters.setting.GaussMea = true;
    parameters.setting.meaIsVec = true;
    parameters.setting.vecRefInertial = false;
    parameters.setting.nVecRef = 1;
    parameters.setting.vRef = [1;0;0];
    parameters.meaNoise = 0.16;
    parameters.initValue.RNoise = eye(3)/200;
    parameters.initValue.xNoise = 0.1^2*eye(3);
    parameters.initValue.U = eye(3);
    parameters.initValue.V = eye(3);
    parameters.initValue.Miu = [0;0;0];

    [gyro,RMea,RTrue,xTrue] = genTrig(t,sf,parameters);
    
    parameters.initValue.U = RTrue(:,:,1)*expRot([0;pi;0]);
    parameters.initValue.Miu = [0.2;0.2;0.2];

    [RMEKF,xMEKF,SigmaMEKF,TMEKF] = MEKF(gyro,RMea,sf,parameters);
    [RMFGAQS,MFGAQS,TMFGAQS] = MFGAnalytic(gyro,RMea,sf,true,parameters);
    [RMFGUQS,MFGUQS,TMFGUQS] = MFGUnscented(gyro,RMea,sf,true,parameters);
%     [RMFGASQ,MFGASQ,TMFGASQ] = MFGAnalytic(gyro,RMea,sf,false,parameters);
%     [RMFGUSQ,MFGUSQ,TMFGUSQ] = MFGUnscented(gyro,RMea,sf,false,parameters);

    parsave(n,gyro,RMea,RTrue,xTrue,parameters,RMEKF,xMEKF,SigmaMEKF,TMEKF,...
        RMFGAQS,MFGAQS,TMFGAQS,RMFGUQS,MFGUQS,TMFGUQS);

end

rmpath('..\Generate-Path');
rmpath('..');
rmpath('..\Matrix-Fisher-Distribution');
rmpath('..\..\rotation3d');

end


function [] = parsave(n,gyro,RMea,RTrue,xTrue,parameters,RMEKF,xMEKF,SigmaMEKF,TMEKF,...
    RMFGAQS,MFGAQS,TMFGAQS,RMFGUQS,MFGUQS,TMFGUQS)

save(strcat('D:\result-SO3Euclid\6-28-2020\',num2str(n)),'gyro','RMea','RTrue','xTrue','parameters','RMEKF','xMEKF','SigmaMEKF','TMEKF',...
        'RMFGAQS','MFGAQS','TMFGAQS','RMFGUQS','MFGUQS','TMFGUQS');

end

