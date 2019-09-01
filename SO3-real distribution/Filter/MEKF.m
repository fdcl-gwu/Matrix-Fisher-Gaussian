function [ R, Sigma ] = MEKF( gyro, RInit, RMea, sf )

filePath = mfilename('fullpath');
pathCell = regexp(path, pathsep, 'split');
if ~any(strcmp(pathCell,getAbsPath('..\..\rotation3d',filePath)))
    addpath(getAbsPath('..\..\rotation3d',filePath));
end

dt = 1/sf;
N = size(gyro,2);

% noise parameters
randomWalk = 10*pi/180;
rotMeaNoise = 0.2;

% data containers
R = zeros(3,3,N);
Sigma = zeros(3,3,N);

% initialize
R(:,:,1) = RInit;
Sigma(:,:,1) = eye(3)*0.01^2;

% filter iteration
for n = 2:N
    % propagate
    av = (gyro(:,n-1)+gyro(:,n))/2;
    F = expRot(av*dt)';
    R(:,:,n) = R(:,:,n-1)*expRot(av*dt);
    Sigma(:,:,n) = F*Sigma(:,:,n-1)*F'+eye(3)*randomWalk^2*dt;
    
    % update
    K = Sigma(:,:,n)*(Sigma(:,:,n)+eye(3)*rotMeaNoise^2)^-1;
    dx = K*logRot(R(:,:,n)'*RMea(:,:,n),'v');
    R(:,:,n) = R(:,:,n)*expRot(dx);
    Sigma(:,:,n) = (eye(3)-K)*Sigma(:,:,n);
end

if ~any(strcmp(pathCell,getAbsPath('..\..\rotation3d',filePath)))
    rmpath(getAbsPath('..\..\rotation3d',filePath));
end

end

