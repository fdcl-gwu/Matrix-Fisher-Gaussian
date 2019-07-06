function [  ] = SO3RealDistThreePara(  )

% Parameters
n1 = 1;
n2 = 3;
p = 3;

Miu1 = 0;
Miu2 = mat2vec(eye(3));
Sigma1 = 1;
K = eye(3)*10;
PReduced = [0,0.7,0]/sqrt(10);  % should be changed

% tangent space
t3 = vec2mat(Miu2)*[0, -1, 0
                    1, 0, 0
                    0, 0, 0]/sqrt(2);
t2 = vec2mat(Miu2)*[0, 0, 1
                    0, 0, 0
                    -1, 0, 0]/sqrt(2);
t1 = vec2mat(Miu2)*[0, 0, 0
                    0, 0, -1
                    0, 1, 0]/sqrt(2);
t1 = mat2vec(t1);
t2 = mat2vec(t2);
t3 = mat2vec(t3);
n = null([t1';t2';t3']);
Rt = [t1';t2';t3';n'];

P = [PReduced,zeros(1,6)]*Rt;

% intermediate parameters
Sigma2 = [K^-1,zeros(3),zeros(3)
          zeros(3),K^-1,zeros(3)
          zeros(3),zeros(3),K^-1];
F = vec2mat(Miu2)*K;

Miuc = @(x2)Miu1+P*Sigma2^-1*(x2-Miu2);
Sigmac = Sigma1-P*Sigma2^-1*P';

% Normalizing constant
c1 = 1/sqrt((2*pi)^n1*det(Sigmac));
c2 = 1;  % should be changed

% density
f = @(x1,x2)1/c1*exp(-1/2*(x1-Miuc(x2))'*Sigmac^-1*(x1-Miuc(x2)))*...
    1/c2*exp(trace(F'*vec2mat(x2)));

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
    fRot = @(x2)f(x1(nx1),mat2vec(x2));
    c(:,:,nx1) = plotRotDist(theta1,theta2,fRot);
end
cmax = max(max(max(c)));

% plot
for nx1 = 1:Nx1
    figure;
    surf(s1,s2,s3,c(:,:,nx1),'LineStyle','none');
    axis equal;
    view([1,1,1]);
    caxis([0,cmax]);
end

end


function [ vec ] = mat2vec( mat )

vec = [mat(1,:)'; mat(2,:)'; mat(3,:)'];

end


function [ mat ] = vec2mat( vec )

mat = [vec(1:3)'; vec(4:6)'; vec(7:9)'];

end

