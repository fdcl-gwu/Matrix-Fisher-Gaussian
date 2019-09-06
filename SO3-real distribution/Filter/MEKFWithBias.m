function [ R, x, G ] = MEKFWithBias( gyro, RInit, RMea, sf )

filePath = mfilename('fullpath');
pathCell = regexp(path, pathsep, 'split');
if ~any(strcmp(pathCell,getAbsPath('..\..\rotation3d',filePath)))
    addpath(getAbsPath('..\..\rotation3d',filePath));
end

dt = 1/sf;
N = size(gyro,2);

% noise parameters
randomWalk = 10*pi/180;
biasInstability = 500/3600*pi/180;
rotMeaNoise = 0.05;

% initialize distribution
Sigma = [eye(3)*100,zeros(3);zeros(3),eye(3)*0.1^2];

% data containers
G.Sigma = zeros(6,6,N); G.Sigma(:,:,1) = Sigma;
R = zeros(3,3,N); R(:,:,1) = RInit;
x = zeros(3,N);

% filter iteration
for n = 2:N
    % propagate
    av = (gyro(:,n-1)+gyro(:,n))/2-x(:,n-1);
    Rp = R(:,:,n-1)*expRot(av*dt);
    F = [expRot(av*dt)',-eye(3)*dt;zeros(3),eye(3)];
    Sigma = F*Sigma*F'+[eye(3)*randomWalk^2*dt,zeros(3);
        zeros(3),eye(3)*biasInstability^2*dt];
    
    % update
    H = [eye(3),zeros(3)];
    K = Sigma*H'*(H*Sigma*H'+eye(3)*rotMeaNoise^2)^-1;
    dx = K*logRot(Rp'*RMea(:,:,n),'v');
    Sigma = (eye(6)-K*H)*Sigma;
    
    R(:,:,n) = Rp*expRot(dx(1:3));
    x(:,n) = x(:,n-1)+dx(4:6);
    
    % record covariance
    G.Sigma(:,:,n) = Sigma;
end

if ~any(strcmp(pathCell,getAbsPath('..\..\rotation3d',filePath)))
    rmpath(getAbsPath('..\..\rotation3d',filePath));
end

end

