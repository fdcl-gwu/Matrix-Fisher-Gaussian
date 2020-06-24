function [] = test_attMea()

addpath('..\..\SO3-real distribution\Generate-Path');
addpath('..');
addpath('..\..\SO3-real distribution\Matrix-Fisher-Distribution');
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
    parameters.t = 60;
    parameters.dt = 1/sf;
    parameters.randomWalk = 10*pi/180;
    parameters.biasInstability = 500/3600*pi/180;
    parameters.GaussMea = false;
    parameters.meaIsVec = false;
    parameters.attMeaLocal = true;
    parameters.rotMeaNoise = diag([2.4,2.4,2.4]);
    parameters.initRNoise = diag([2.4,2.4,2.4]);
    parameters.initXNoise = 0.1^2*eye(3);
    parameters.RInit = eye(3);
    parameters.xInit = [0;0;0];

    [gyro,RMea,RTrue,xTrue] = genTrig_attMea(t,sf,parameters);

    % stachastic settings
    parameters.RInit = RMea(:,:,1);
    parameters.xInit = [0;0;0];

    [RMEKF,xMEKF,G,TMEKF] = MEKF(gyro,RMea,parameters);
    [RMFGA,MFGA,TMFGA] = MFGAnalytic(gyro,RMea,parameters);
    [RMFGU,MFGU,TMFGU] = MFGUnscented(gyro,RMea,parameters);

    parsave(n,gyro,RMea,RTrue,xTrue,parameters,RMEKF,xMEKF,G,TMEKF,...
        RMFGA,MFGA,TMFGA,RMFGU,MFGU,TMFGU);

end

rmpath('..\..\SO3-real distribution\Generate-Path');
rmpath('..');
rmpath('..\..\SO3-real distribution\Matrix-Fisher-Distribution');
rmpath('..\..\rotation3d');

end


function [] = parsave(n,gyro,RMea,RTrue,xTrue,parameters,RMEKF,xMEKF,G,TMEKF,...
    RMFGA,MFGA,TMFGA,RMFGU,MFGU,TMFGU)

save(strcat('D:\result-SO3Euclid\6-21-2020\',num2str(n)),'gyro','RMea','RTrue','xTrue','parameters','RMEKF','xMEKF','G','TMEKF',...
        'RMFGA','MFGA','TMFGA','RMFGU','MFGU','TMFGU');

end

