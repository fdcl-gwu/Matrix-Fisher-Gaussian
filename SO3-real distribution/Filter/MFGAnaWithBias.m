function [ R, MFG ] = MFGAnaWithBias( gyro, RInit, RMea, sf )

filePath = mfilename('fullpath');
pathCell = regexp(path, pathsep, 'split');
if ~any(strcmp(pathCell,getAbsPath('..\..\rotation3d',filePath)))
    addpath(getAbsPath('..\..\rotation3d',filePath));
end
if ~any(strcmp(pathCell,getAbsPath('..\Matrix-Fisher-Distribution',filePath)))
    addpath(getAbsPath('..\Matrix-Fisher-Distribution',filePath));
end
if ~any(strcmp(pathCell,getAbsPath('..\',filePath)))
    addpath(getAbsPath('..\',filePath));
end

dt = 1/sf;
N = size(gyro,2);

% noise parameters
randomWalk = 10*pi/180;
biasInstability = 500/3600*pi/180;
rotMeaNoise = 0.2;

SM = Gau2MF(rotMeaNoise);

% initialize distribution
Miu = [0;0;0];
Sigma = eye(3)*0.1^2;
P = zeros(3);
U = RInit*expRot([pi,0,0]);
V = eye(3);
S = Gau2MF(sqrt(1/200));
S(2,2) = S(2,2)+1e-5;
S(3,3) = S(3,3)+2e-5;

% data containers
MFG.Miu = zeros(3,N); MFG.Miu(:,1) = Miu;
MFG.Sigma = zeros(3,3,N); MFG.Sigma(:,:,1) = Sigma;
MFG.P = zeros(3,3,N); MFG.P(:,:,1) = P;
MFG.U = zeros(3,3,N); MFG.U(:,:,1) = U;
MFG.V = zeros(3,3,N); MFG.V(:,:,1) = V;
MFG.S = zeros(3,N); MFG.S(:,1) = diag(S);
R = zeros(3,3,N); R(:,:,1) = RInit*expRot([pi,0,0]);

% filter iteration
for n = 2:N
    % uncertainty propagation
    omega = (gyro(:,n-1)+gyro(:,n))/2;
    [Miu,Sigma,P,U,S,V] = MFGGyroProp(omega,Miu,Sigma,P,U,S,V,randomWalk*eye(3),biasInstability^2*dt*eye(3),dt);
    
    % update
    if rem(n,5)==0
        [Miu,Sigma,P,U,S,V] = MFGMulMF(Miu,Sigma,P,U,S,V,RMea(:,:,n)*SM);
    end
    
    % record results
    MFG.Miu(:,n) = Miu;
    MFG.Sigma(:,:,n) = Sigma;
    MFG.P(:,:,n) = P;
    MFG.U(:,:,n) = U;
    MFG.V(:,:,n) = V;
    MFG.S(:,n) = diag(S);
    R(:,:,n) = U*V';
end

if ~any(strcmp(pathCell,getAbsPath('..\..\rotation3d',filePath)))
    rmpath(getAbsPath('..\..\rotation3d',filePath));
end
if ~any(strcmp(pathCell,getAbsPath('..\Matrix-Fisher-Distribution',filePath)))
    addpath(getAbsPath('..\Matrix-Fisher-Distribution',filePath));
end
if ~any(strcmp(pathCell,getAbsPath('..\',filePath)))
    rmpath(getAbsPath('..\',filePath));
end

end


function [ S ] = Gau2MF( sigma )

N = 100000;
v = sigma*randn(3,N);

R = expRot(v);
ER = mean(R,3);
[~,D,~] = psvd(ER);

S = pdf_MF_M2S(diag(D));
S = eye(3)*mean(S);

end

