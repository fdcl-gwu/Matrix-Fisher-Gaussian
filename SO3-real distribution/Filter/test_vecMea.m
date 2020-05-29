function [  ] = test_vecMea(  )

addpath('..\Generate-Path');
addpath('..');
addpath('..\Matrix-Fisher-Distribution');
addpath('..\..\rotation3d');

rng(1);

t = 600;
sf = 150;


% general settings
parameters = [];
parameters.t = t;
parameters.dt = 1/sf;
parameters.randomWalk = 10*pi/180;
parameters.biasInstability = 500/3600*pi/180;
parameters.GaussMea = false;
parameters.meaIsVec = true;
parameters.meaNoise = 200;
parameters.initRNoise = diag([200,0,0]);
parameters.initXNoise = 0.1^2*eye(3);
parameters.RInit = eye(3);
parameters.xInit = [0;0;0];

[gyro,RTrue,xTrue,vMea,vRef] = genTrig_vecMea(t,sf,parameters);

% stachastic settings
parameters.RInit = eye(3);
parameters.xInit = [0;0;0];

[RMEKF,xMEKF,G,TMEKF] = MEKF(gyro,{vMea,vRef},parameters);
[RMFGA,MFGA,TMFGA] = MFGAnalytic(gyro,{vMea,vRef},parameters);
[RMFGU,MFGU,TMFGU] = MFGUnscented(gyro,{vMea,vRef},parameters);

save('D:\result-SO3Euclid\5-16-2020','gyro','vRef','vMea','RTrue','xTrue','parameters','RMEKF','xMEKF','G','TMEKF','RMFGA','MFGA','TMFGA','RMFGU','MFGU','TMFGU');

rmpath('..\Generate-Path');
rmpath('..');
rmpath('..\Matrix-Fisher-Distribution');
rmpath('..\..\rotation3d');

end

