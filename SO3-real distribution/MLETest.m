function [  ] = MLETest(  )
addpath('Matrix-Fisher-Distribution');
addpath('..\rotation3d');

% parameters
N = 1;
Ns = 100000;

Miu = 0;
Sigma = 1;
U = expRM([0,0,0]);
V = expRM([-0.3,0.9,0.5]);
S = diag([25,25,-25]);
F = U*S*V';
PReduced = [0.5,0.5,0.5]/sqrt(25);

% intermediate parameters
[U,S,V] = usvd(F,true);
M = U*V';
K = V*S*V';
Miu2 = mat2vec(M);
Sigma2Inv = [K,zeros(3),zeros(3)
             zeros(3),K,zeros(3)
             zeros(3),zeros(3),K];

% tangent space
[U0,~,V0] = usvd(F,false);
M0 = U0*V0';
Omega1 = skew(V0*[1;0;0]);
Omega2 = skew(V0*[0;1;0]);
Omega3 = skew(V0*[0;0;1]);
t1 = mat2vec(M0*Omega1)/sqrt(2);
t2 = mat2vec(M0*Omega2)/sqrt(2);
t3 = mat2vec(M0*Omega3)/sqrt(2);
n = null([t1';t2';t3']);
Rt = [t1';t2';t3';n'];

P = [PReduced,zeros(1,6)]*Rt;
         
% sample Matrix Fisher distribution
Rs = pdf_MF_sampling(F,Ns);

% sample from conditional Gaussian distribution
xs = zeros(N,Ns);
Sigmac = Sigma-P*Sigma2Inv*P';
for ns = 1:Ns
    Miuc = Miu+P*Sigma2Inv*(mat2vec(Rs(:,:,ns))-Miu2);
    xs(:,ns) = mvnrnd(Miuc,Sigmac,1);
end

% maximum likelihood estimation
[MiuHat,SigmaHat,PHat,MHat,KHat,PTildeHat] = SO3RealMLE(xs,Rs);

rmpath('Matrix-Fisher-Distribution');
rmpath('..\rotation3d');

end


function [ vec ] = mat2vec( mat )

vec = [mat(1,:)'; mat(2,:)'; mat(3,:)'];

end

