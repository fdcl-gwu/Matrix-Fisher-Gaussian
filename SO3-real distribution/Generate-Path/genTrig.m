function [ gyroMea, RMea, RTrue, biasTrue ] = genTrig( t, sf, parameters )

filePath = mfilename('fullpath');
pathCell = regexp(path, pathsep, 'split');
if ~any(strcmp(pathCell,getAbsPath('..\..\rotation3d',filePath)))
    addpath(getAbsPath('..\..\rotation3d',filePath));
end
if ~any(strcmp(pathCell,getAbsPath('..\Matrix-Fisher-Distribution',filePath)))
    addpath(getAbsPath('..\Matrix-Fisher-Distribution',filePath));
end

time = (0:1/sf:t);
N = length(time);

% motion parameters
E.fr = 0.35; E.fp = 0.35; E.fh = 0.35;
E.magr = pi; E.magp = pi/2; E.magh = pi;

% noise parameters
if exist('parameters','var')
    randomWalk = parameters.randomWalk;
    biasInstability = parameters.biasInstability;
    rotMeaNoise = parameters.rotMeaNoise;
else
    randomWalk = 10*pi/180;
    biasInstability = 500/3600*pi/180;
    rotMeaNoise = 0.2^2*eye(3);
end

% true state
roll = @(t)E.magr*sin(E.fr*2*pi*t);
pitch = @(t)E.magp*sin(E.fp*2*pi*t);
yaw = @(t)E.magh*sin(E.fh*2*pi*t);

sr2 = @(t)sin(roll(t)/2); sp2 = @(t)sin(pitch(t)/2); sy2 = @(t)sin(yaw(t)/2);
cr2 = @(t)cos(roll(t)/2); cp2 = @(t)cos(pitch(t)/2); cy2 = @(t)cos(yaw(t)/2);
q1 = @(t)cr2(t).*cp2(t).*cy2(t)+sr2(t).*sp2(t).*sy2(t);
q2 = @(t)sr2(t).*cp2(t).*cy2(t)-cr2(t).*sp2(t).*sy2(t);
q3 = @(t)cr2(t).*sp2(t).*cy2(t)+sr2(t).*cp2(t).*sy2(t);
q4 = @(t)cr2(t).*cp2(t).*sy2(t)-sr2(t).*sp2(t).*cy2(t);
qTrue = [q1(time);q2(time);q3(time);q4(time)];
RTrue = qua2rot(qTrue);

dr = @(t)E.magr*E.fr*2*pi*cos(E.fr*2*pi*t);
dp = @(t)E.magp*E.fp*2*pi*cos(E.fp*2*pi*t);
dy = @(t)E.magh*E.fh*2*pi*cos(E.fh*2*pi*t);

dq1 = @(t)-sr2(t).*cp2(t).*cy2(t)*0.5.*dr(t)-cr2(t).*sp2(t).*cy2(t)*0.5.*dp(t)-cr2(t).*cp2(t).*sy2(t)*0.5.*dy(t)...
          +cr2(t).*sp2(t).*sy2(t)*0.5.*dr(t)+sr2(t).*cp2(t).*sy2(t)*0.5.*dp(t)+sr2(t).*sp2(t).*cy2(t)*0.5.*dy(t);
dq2 = @(t)+cr2(t).*cp2(t).*cy2(t)*0.5.*dr(t)-sr2(t).*sp2(t).*cy2(t)*0.5.*dp(t)-sr2(t).*cp2(t).*sy2(t)*0.5.*dy(t)...
          +sr2(t).*sp2(t).*sy2(t)*0.5.*dr(t)-cr2(t).*cp2(t).*sy2(t)*0.5.*dp(t)-cr2(t).*sp2(t).*cy2(t)*0.5.*dy(t);
dq3 = @(t)-sr2(t).*sp2(t).*cy2(t)*0.5.*dr(t)+cr2(t).*cp2(t).*cy2(t)*0.5.*dp(t)-cr2(t).*sp2(t).*sy2(t)*0.5.*dy(t)...
          +cr2(t).*cp2(t).*sy2(t)*0.5.*dr(t)-sr2(t).*sp2(t).*sy2(t)*0.5.*dp(t)+sr2(t).*cp2(t).*cy2(t)*0.5.*dy(t);
dq4 = @(t)-sr2(t).*cp2(t).*sy2(t)*0.5.*dr(t)-cr2(t).*sp2(t).*sy2(t)*0.5.*dp(t)+cr2(t).*cp2(t).*cy2(t)*0.5.*dy(t)...
          -cr2(t).*sp2(t).*cy2(t)*0.5.*dr(t)-sr2(t).*cp2(t).*cy2(t)*0.5.*dp(t)+sr2(t).*sp2(t).*sy2(t)*0.5.*dy(t);

w1 = @(t)(q1(t).*dq1(t)+q2(t).*dq2(t)+q3(t).*dq3(t)+q4(t).*dq4(t))*2;
w2 = @(t)(q1(t).*dq2(t)-q2(t).*dq1(t)-q3(t).*dq4(t)+q4(t).*dq3(t))*2;
w3 = @(t)(q1(t).*dq3(t)+q2(t).*dq4(t)-q3(t).*dq1(t)-q4(t).*dq2(t))*2;
w4 = @(t)(q1(t).*dq4(t)-q2(t).*dq3(t)+q3(t).*dq2(t)-q4(t).*dq1(t))*2;

gyro = [w2(time);w3(time);w4(time)];

% add noise
biasNoise = randn(3,N)*biasInstability*sqrt(sf);
biasTrue = cumsum(biasNoise/sf,2);

gyroNoise = randn(3,N)*randomWalk*sqrt(sf);
gyroMea = gyro+gyroNoise+biasTrue;

if parameters.GaussMea
    RNoise = expRot(mvnrnd([0;0;0],rotMeaNoise,N));
else
    RNoise = pdf_MF_sampling(rotMeaNoise,N);
end
RMea = zeros(3,3,N);
for n = 1:N
    RMea(:,:,n) = RTrue(:,:,n)*RNoise(:,:,n)';
end

if ~any(strcmp(pathCell,getAbsPath('..\..\rotation3d',filePath)))
    rmpath(getAbsPath('..\..\rotation3d',filePath));
end
if ~any(strcmp(pathCell,getAbsPath('..\Matrix-Fisher-Distribution',filePath)))
    rmpath(getAbsPath('..\Matrix-Fisher-Distribution',filePath));
end

end

