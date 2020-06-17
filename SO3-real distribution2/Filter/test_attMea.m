function [] = test_attMea()

addpath('..\Generate-Path');
addpath('..');
addpath('..\Matrix-Fisher-Distribution');
addpath('..\..\rotation3d');

t = 600;
sf = 150;

parameters = [];
parameters.t = 60;
parameters.dt = 1/sf;
parameters.randomWalk = 10*pi/180;
parameters.biasInstability = 500/3600*pi/180;
parameters.GaussMea = false;
parameters.meaIsVec = false;
parameters.attMeaLocal = true;
parameters.rotMeaNoise = diag([200,0,0]);
parameters.initRNoise = diag([200,200,200]);
parameters.initXNoise = 0.1^2*eye(3);
parameters.RInit = eye(3);
parameters.xInit = [0;0;0];

[gyro,RMea,RTrue,xTrue] = genTrig_attMea(t,sf,parameters);

% stachastic settings
parameters.RInit = eye(3);
parameters.xInit = [0;0;0];

[RMEKF,xMEKF,G,TMEKF] = MEKF(gyro,RMea,parameters);
[RMFGA,MFGA,TMFGA] = MFGAnalytic(gyro,RMea,parameters);
[RMFGU,MFGU,TMFGU] = MFGUnscented(gyro,RMea,parameters);

save('D:\result-SO3Euclid\5-15-2020','gyro','RMea','RTrue','xTrue','parameters','RMEKF','xMEKF','G','TMEKF','RMFGA','MFGA','TMFGA','RMFGU','MFGU','TMFGU');

rmpath('..\Generate-Path');
rmpath('..');
rmpath('..\Matrix-Fisher-Distribution');
rmpath('..\..\rotation3d');

end

