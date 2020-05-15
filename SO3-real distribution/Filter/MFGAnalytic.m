function [ R, MFG, stepT ] = MFGAnalytic( gyro, mea, parameters )

N = size(gyro,2);
dt = parameters.dt;

% noise parameters
randomWalk = parameters.randomWalk;
biasInstability = parameters.biasInstability;
if parameters.meaIsVec
    vecMeaNoise = parameters.meaNoise;
else
    if parameters.GaussMea
        SM = Gau2MF(parameters.rotMeaNoise);
    else
        SM = parameters.rotMeaNoise;
    end
end

% measurement
if parameters.meaIsVec
    vMea = mea{1};
    vRef = mea{2};
end

% initialize distribution
Miu = -parameters.xInit;
Sigma = parameters.initXNoise;
P = zeros(3);
U = parameters.RInit;
V = eye(3);
if parameters.GaussMea
    S = Gau2MF(parameters.initRNoise);
else
    S = parameters.initRNoise;
end
S(2,2) = S(2,2)+1e-5;
S(3,3) = S(3,3)+2e-5;

% data containers
MFG.Miu = zeros(3,N); MFG.Miu(:,1) = Miu;
MFG.Sigma = zeros(3,3,N); MFG.Sigma(:,:,1) = Sigma;
MFG.P = zeros(3,3,N); MFG.P(:,:,1) = P;
MFG.U = zeros(3,3,N); MFG.U(:,:,1) = U;
MFG.V = zeros(3,3,N); MFG.V(:,:,1) = V;
MFG.S = zeros(3,N); MFG.S(:,1) = diag(S);
R = zeros(3,3,N); R(:,:,1) = U*V';
stepT = zeros(N-1,1);

% filter iteration
try
for n = 2:N
    tic;
    % uncertainty propagation
    omega = (gyro(:,n-1)+gyro(:,n))/2;
    [Miu,Sigma,P,U,S,V] = MFGGyroProp(omega,Miu,Sigma,P,U,S,V,randomWalk*eye(3),biasInstability^2*dt*eye(3),dt);
    
    % update
    if rem(n,5)==0
        if parameters.meaIsVec
            if size(vRef,1)==3
                [Miu,Sigma,P,U,S,V] = MFGMulMF(Miu,Sigma,P,U,S,V,vecMeaNoise*(vRef(:,n)*vMea(:,n)'));
            else
                [Miu,Sigma,P,U,S,V] = MFGMulMF(Miu,Sigma,P,U,S,V,vecMeaNoise*(vRef(1:3,n)*vMea(1:3,n)'+vRef(4:6,n)*vMea(4:6,n)'));
            end
        else
            [Miu,Sigma,P,U,S,V] = MFGMulMF(Miu,Sigma,P,U,S,V,mea(:,:,n)*SM);
        end
    end
    
    % record results
    MFG.Miu(:,n) = Miu;
    MFG.Sigma(:,:,n) = Sigma;
    MFG.P(:,:,n) = P;
    MFG.U(:,:,n) = U;
    MFG.V(:,:,n) = V;
    MFG.S(:,n) = diag(S);
    R(:,:,n) = U*V';
    
    stepT(n-1) = toc;
end
catch
    pause(1);
end

end


function [ S ] = Gau2MF( Sigma )

N = 100000;
v = mvnrnd([0;0;0],Sigma,N);

R = expRot(v);
ER = mean(R,3);
[~,D,~] = psvd(ER);

S = diag(pdf_MF_M2S(diag(D)));

end

