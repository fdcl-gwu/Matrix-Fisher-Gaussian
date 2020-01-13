function [ R, x, G, stepT ] = MEKFWithBias( gyro, RMea, parameters )

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
rotMeaNoise = parameters.rotMeaNoise;

% initialize distribution
Sigma = [eye(3)*parameters.initRsigma^2,zeros(3)
    zeros(3),eye(3)*parameters.initXsigma^2];

% data containers
G.Sigma = zeros(6,6,N); G.Sigma(:,:,1) = Sigma;
R = zeros(3,3,N); R(:,:,1) = parameters.RInit;
x = zeros(3,N); x(:,1) = parameters.xInit;
stepT = zeros(N-1,1);

% filter iteration
for n = 2:N
    tic;
    
    % propagate
    av = (gyro(:,n-1)+gyro(:,n))/2-x(:,n-1);
    Rp = R(:,:,n-1)*expRot(av*dt);
    F = [expRot(av*dt)',-eye(3)*dt;zeros(3),eye(3)];
    Sigma = F*Sigma*F'+[eye(3)*randomWalk^2*dt,zeros(3);
        zeros(3),eye(3)*biasInstability^2*dt];
    
    % update
    if rem(n,5)==0
        H = [eye(3),zeros(3)];
        K = Sigma*H'*(H*Sigma*H'+eye(3)*rotMeaNoise^2)^-1;
        dx = K*logRot(Rp'*RMea(:,:,n),'v');
        Sigma = (eye(6)-K*H)*Sigma;

        R(:,:,n) = Rp*expRot(dx(1:3));
        x(:,n) = x(:,n-1)+dx(4:6);
    else
        R(:,:,n) = Rp;
        x(:,n) = x(:,n-1);
    end
    
    % record covariance
    G.Sigma(:,:,n) = Sigma;
    
    stepT(n-1) = toc;
end

if ~any(strcmp(pathCell,getAbsPath('..\..\rotation3d',filePath)))
    rmpath(getAbsPath('..\..\rotation3d',filePath));
end

end

