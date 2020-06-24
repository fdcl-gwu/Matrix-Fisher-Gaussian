function [ R, x, MFG, stepT ] = MFGAnalyticIMU( gyro, acce, RMea, pMea, parameters, RTrue, xTrue )

N = size(gyro,2);
dt = parameters.dt;

% noise parameters
randomWalk = parameters.randomWalk;
biasInstability = parameters.biasInstability;
acceRandomWalk = parameters.acceRandomWalk;
acceBiasInstability = parameters.acceBiasInstability;
if parameters.GaussMea
    SM = Gau2MF(parameters.rotMeaNoise);
else
    SM = parameters.rotMeaNoise;
end
posMeaNoise = parameters.posMeaNoise;

H.Hgu = randomWalk*eye(3);
H.Hgv = biasInstability*eye(3);
H.Hau = acceRandomWalk*eye(3);
H.Hav = acceBiasInstability*eye(3);

% initialize distribution
Miu = parameters.xInit;
Sigma = parameters.initXNoise;
P = zeros(12,3);
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
MFG.Miu = zeros(12,N); MFG.Miu(:,1) = Miu;
MFG.Sigma = zeros(12,12,N); MFG.Sigma(:,:,1) = Sigma;
MFG.P = zeros(12,3,N); MFG.P(:,:,1) = P;
MFG.U = zeros(3,3,N); MFG.U(:,:,1) = U;
MFG.V = zeros(3,3,N); MFG.V(:,:,1) = V;
MFG.S = zeros(3,N); MFG.S(:,1) = diag(S);
R = zeros(3,3,N); R(:,:,1) = U*V';
x = zeros(12,N); x(:,1) = Miu;
stepT = zeros(N-1,1);

% filter iteration
for n = 2:N
    tic;
    
    % uncertainty propagation
    omega = (gyro(:,n-1)+gyro(:,n))/2;
    a = (acce(:,n-1)+acce(:,n))/2;
    [Miu,Sigma,P,U,S,V] = MFGIMUProp(omega,a,Miu,Sigma,P,U,S,V,H,dt);
    
    % unscented update
    if rem(n,5)==0
%         [xs,Rs] = MFGSampling(Miu,Sigma,P,U,S,V,10000);
%         w = ones(1,10000)/10000;
% 
%         for i = 1:10000
%             dR = Rs(:,:,i)'*RMea(:,:,n);
%             dx = xs(4:6,i)-pMea(:,n);
%             w(i) = w(i)*exp(trace(SM*dR))*exp(-0.5*dx'*posMeaNoise^-1*dx);
%         end
%         
%         w = w/sum(w);
%         [Miu,Sigma,P,U,S,V] = MFGMLEAppro(xs,Rs,w);
        
    end
    
    % record results
    MFG.Miu(:,n) = Miu;
    MFG.Sigma(:,:,n) = Sigma;
    MFG.P(:,:,n) = P;
    MFG.U(:,:,n) = U;
    MFG.V(:,:,n) = V;
    MFG.S(:,n) = diag(S);
    R(:,:,n) = U*V';
    x(:,n) = Miu;
    
    stepT(n-1) = toc;
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
