function [ R, x, G, stepT ] = MEKFIMU( gyro, acce, RMea, pMea, parameters)

filePath = mfilename('fullpath');
pathCell = regexp(path, pathsep, 'split');
if ~any(strcmp(pathCell,getAbsPath('..\..\rotation3d',filePath)))
    addpath(getAbsPath('..\..\rotation3d',filePath));
end

N = size(gyro,2);
dt = parameters.dt;

% noise parameters
randomWalk = parameters.randomWalk;
biasInstability = parameters.biasInstability;
acceRandomWalk = parameters.acceRandomWalk;
acceBiasInstability = parameters.acceBiasInstability;
if parameters.GaussMea
    rotMeaNoise = parameters.rotMeaNoise;
else
    rotMeaNoise = MF2Gau(parameters.rotMeaNoise);
end
posMeaNoise = parameters.posMeaNoise;

% initialize distribution
if parameters.GaussMea
    initRNoise = parameters.initRNoise;
else
    initRNoise = MF2Gau(parameters.initRNoise);
end
Sigma = [initRNoise, zeros(3,12);
    zeros(12,3), parameters.initXNoise];

% filter parameters
Q = [eye(3)*randomWalk^2*dt,zeros(3,12);
    zeros(3),eye(3)*biasInstability^2*dt,zeros(3,9);
    zeros(3,15);
    zeros(3,9),eye(3)*acceRandomWalk^2*dt,zeros(3);
    zeros(3,12),eye(3)*acceBiasInstability^2*dt];
meaNoise = [rotMeaNoise,zeros(3);
    zeros(3),posMeaNoise];
H = [eye(3),zeros(3,12);
    zeros(3,6),eye(3),zeros(3,6)];

% data containers
G.Sigma = zeros(15,15,N); G.Sigma(:,:,1) = Sigma;
R = zeros(3,3,N); R(:,:,1) = parameters.RInit;
x = zeros(12,N); x(:,1) = parameters.xInit;
stepT = zeros(N-1,1);

% filter iteration
for n = 2:N
    tic;
    
    % propagate
    av = (gyro(:,n-1)+gyro(:,n))/2-x(1:3,n-1);
    a = (acce(:,n-1)+acce(:,n))/2-x(10:12,n-1);
    Rp = R(:,:,n-1)*expRot(av*dt);
    pp = x(4:6,n-1) + dt*x(7:9,n-1);
    vp = x(7:9,n-1) + dt*(R(:,:,n-1)*a-[0;0;9.8]);
    xp = [x(1:3,n-1);pp;vp;x(10:12,n-1)];
    
    F = [expRot(av*dt)',-eye(3)*dt,zeros(3,9);
        zeros(3,3),eye(3),zeros(3,9);
        zeros(3,6),eye(3),eye(3)*dt,zeros(3,3);
        -dt*R(:,:,n-1)*hat(a),zeros(3,6),eye(3),-R(:,:,n-1)*dt;
        zeros(3,12),eye(3)];
    Sigma = F*Sigma*F'+Q;
    
    % update
    if rem(n,5)==0
        K = Sigma*H'*(H*Sigma*H'+meaNoise)^-1;
        y = [logRot(Rp'*RMea(:,:,n),'v');pMea(:,n)-xp(4:6)];
        dx = K*y;
        Sigma = (eye(15)-K*H)*Sigma;

        R(:,:,n) = Rp*expRot(dx(1:3));
        x(:,n) = xp+dx(4:15);
    else
        R(:,:,n) = Rp;
        x(:,n) = xp;
    end
    
    % record covariance
    G.Sigma(:,:,n) = Sigma;
    
    stepT(n-1) = toc;
end

end


function [ Sigma ] = MF2Gau( S )

filePath = mfilename('fullpath');
pathCell = regexp(path, pathsep, 'split');
if ~any(strcmp(pathCell,getAbsPath('..\Matrix-Fisher-Distribution',filePath)))
    addpath(getAbsPath('..\Matrix-Fisher-Distribution',filePath));
end

N = 100000;
R = pdf_MF_sampling(S,N);

v = logRot(R,'v');
Sigma = cov(v');

if ~any(strcmp(pathCell,getAbsPath('..\Matrix-Fisher-Distribution',filePath)))
    addpath(getAbsPath('..\Matrix-Fisher-Distribution',filePath));
end

end

