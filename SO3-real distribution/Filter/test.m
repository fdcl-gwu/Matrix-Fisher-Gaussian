function [] = test()

addpath('..\Generate-Path');
addpath('..');
addpath('..\Matrix-Fisher-Distribution');
addpath('..\..\rotation3d');

t = 120;
sf = 150;
N = 100;

% random seeds
parpool(10);
parfor n = 1:10
    rng(n);
end
    
parfor n = 1:N
    parameters = [];
    parameters.setting.omegaLocal = false;
    parameters.setting.gyroFail = false;
    parameters.setting.GaussMea = false;
    parameters.setting.meaIsVec = true;
    parameters.setting.vecRefInertial = true;
    parameters.setting.attMeaLocal = false;
    parameters.setting.nVecRef = 2;
    parameters.setting.vRef = [1;0;0;0;1;0];
    parameters.meaNoise = [1;100];
    parameters.initValue.RNoise = diag([0,0,0]);
    parameters.initValue.xNoise = 0.1^2*eye(3);
    parameters.initValue.U = expRot([pi,0,0]);
    parameters.initValue.V = eye(3);
    parameters.initValue.Miu = [0.2;0.2;0.2];

    [gyro,Mea,RTrue,xTrue] = genTrig(t,sf,parameters);

    try
        [RMEKF,xMEKF,SigmaMEKF,TMEKF] = MEKF(gyro,Mea,sf,parameters);
        [RUKF,xUKF,SigmaUKF,TUKF] = UKF(gyro,Mea,sf,parameters);
        [RMFGAQS,MFGAQS,TMFGAQS] = MFGAnalytic(gyro,Mea,sf,true,parameters);
        [RMFGUQS,MFGUQS,TMFGUQS] = MFGUnscented(gyro,Mea,sf,true,parameters);
        [RMFGASQ,MFGASQ,TMFGASQ] = MFGAnalytic(gyro,Mea,sf,false,parameters);
        [RMFGUSQ,MFGUSQ,TMFGUSQ] = MFGUnscented(gyro,Mea,sf,false,parameters);
    catch
        s = rng;
        disp(strcat('current random seed: ', num2str(s.Seed)));
    end

    parsave(n,gyro,Mea,RTrue,xTrue,parameters,...
        RMEKF,xMEKF,SigmaMEKF,TMEKF,RUKF,xUKF,SigmaUKF,TUKF,...
        RMFGAQS,MFGAQS,TMFGAQS,RMFGUQS,MFGUQS,TMFGUQS,...
        RMFGASQ,MFGASQ,TMFGASQ,RMFGUSQ,MFGUSQ,TMFGUSQ);
end

rmpath('..\Generate-Path');
rmpath('..');
rmpath('..\Matrix-Fisher-Distribution');
rmpath('..\..\rotation3d');

end


function [] = parsave(n,gyro,RMea,RTrue,xTrue,parameters,...
    RMEKF,xMEKF,SigmaMEKF,TMEKF,RUKF,xUKF,SigmaUKF,TUKF,...
    RMFGAQS,MFGAQS,TMFGAQS,RMFGUQS,MFGUQS,TMFGUQS,...
    RMFGASQ,MFGASQ,TMFGASQ,RMFGUSQ,MFGUSQ,TMFGUSQ)

save(strcat('D:\result-SO3Euclid\11-26-2020\',num2str(n)),'gyro','RMea','RTrue','xTrue','parameters',...
        'RMEKF','xMEKF','SigmaMEKF','TMEKF','RUKF','xUKF','SigmaUKF','TUKF',...
        'RMFGAQS','MFGAQS','TMFGAQS','RMFGUQS','MFGUQS','TMFGUQS',...
        'RMFGASQ','MFGASQ','TMFGASQ','RMFGUSQ','MFGUSQ','TMFGUSQ');

end

