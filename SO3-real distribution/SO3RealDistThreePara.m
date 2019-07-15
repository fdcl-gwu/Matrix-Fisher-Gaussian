function [  ] = SO3RealDistThreePara(  )
addpath('Matrix-Fisher-Distribution');
addpath('..\rotation3d');

% parameters
n1 = 1;

Miu = 0;
Sigma = 1;
U = expRM([0,0,0]);
V = expRM([-0.3,0.9,0.5]);
S = diag([25,5,-5]);
PReduced = [0.7,0,0]/sqrt(25);

% intermediate parameters
M = U*V';
K = V*S*V';
F = M*K;
Miu2 = mat2vec(M);
Sigma2Inv = [K,zeros(3),zeros(3)
             zeros(3),K,zeros(3)
             zeros(3),zeros(3),K];

% tangent space
Omega3 = [0, -1, 0
          1, 0, 0
          0, 0, 0]/sqrt(2);
Omega2 = [0, 0, 1
          0, 0, 0
          -1, 0, 0]/sqrt(2);
Omega1 = [0, 0, 0
          0, 0, -1
          0, 1, 0]/sqrt(2);
t1 = mat2vec(M*Omega1);
t2 = mat2vec(M*Omega2);
t3 = mat2vec(M*Omega3);
n = null([t1';t2';t3']);
Rt = [t1';t2';t3';n'];

P = [PReduced,zeros(1,6)]*Rt;

% other intermediate parameters
Miuc = @(R)Miu+P*Sigma2Inv*(mat2vec(R)-Miu2);
Sigmac = Sigma-P*Sigma2Inv*P';

% Normalizing constant
c1 = 1/sqrt((2*pi)^n1*det(Sigmac));
c2 = pdf_MF_normal([S(1,1),S(2,2),S(3,3)]);

% density
f = @(x1,R)1/c1*exp(-1/2*(x1-Miuc(R))'*Sigmac^-1*(x1-Miuc(R)))*...
    1/c2*exp(trace(F*R'));

% spherical grid
Nt1 = 100;
Nt2 = 100;
theta1 = linspace(-pi,pi,Nt2);
theta2 = linspace(0,pi,Nt1);
s1 = cos(theta1)'.*sin(theta2);
s2 = sin(theta1)'.*sin(theta2);
s3 = repmat(cos(theta2),Nt1,1);

% linear grid
Nx1 = 11;
x1 = linspace(-1,1,Nx1);

% color map
c = zeros(Nt1,Nt2,Nx1);
for nx1 = 1:Nx1
    fRot = @(R)f(x1(nx1),R);
    c(:,:,nx1) = plotSO3Dist(theta1,theta2,fRot);
end
cmax = max(max(max(c)));

% plot
for nx1 = 1:Nx1
    figure; hold on;
    surf(s1,s2,s3,c(:,:,nx1),'LineStyle','none');
    plot3([0,M(1,1)*1.1],[0,M(2,1)*1.1],[0,M(3,1)*1.1],'b');
    plot3([0,M(1,2)*1.1],[0,M(2,2)*1.1],[0,M(3,2)*1.1],'r');
    plot3([0,M(1,3)*1.1],[0,M(2,3)*1.1],[0,M(3,3)*1.1],'y');
    
    plot3([0,U(1,1)*1.1],[0,U(2,1)*1.1],[0,U(3,1)*1.1],'b');
    plot3([0,U(1,2)*1.1],[0,U(2,2)*1.1],[0,U(3,2)*1.1],'r');
    plot3([0,U(1,3)*1.1],[0,U(2,3)*1.1],[0,U(3,3)*1.1],'y');
    scatter3(U(1,1)*1.1,U(2,1)*1.1,U(3,1)*1.1,[],'b');
    scatter3(U(1,2)*1.1,U(2,2)*1.1,U(3,2)*1.1,[],'r');
    scatter3(U(1,3)*1.1,U(2,3)*1.1,U(3,3)*1.1,[],'y');
    axis equal;
    view([1,1,1]);
    caxis([0,cmax]);
end

rmpath('Matrix-Fisher-Distribution');
rmpath('..\rotation3d');

end


function [ vec ] = mat2vec( mat )

vec = [mat(1,:)'; mat(2,:)'; mat(3,:)'];

end


function [ mat ] = vec2mat( vec )

mat = [vec(1:3)'; vec(4:6)'; vec(7:9)'];

end

