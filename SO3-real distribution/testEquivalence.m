function [ FError, SigmacError, maxMiucError ] = testEquivalence(  )

pathCell = regexp(path, pathsep, 'split');
if ~any(strcmp(pathCell,getAbsPath('Matrix-Fisher-Distribution')))
    addpath('Matrix-Fisher-Distribution');
end
if ~any(strcmp(pathCell,getAbsPath('..\rotation3d')))
    addpath('..\rotation3d');
end

%% distribution 1
% parameters
Miu1 = 0;
Sigma1 = 1;
U1 = expRM([0.2,0.3,0.4]);
V1 = expRM([-0.1,0.5,0.7]);
S1 = diag([25,10,-10]);
PTilde1 = [0.5,-0.4,0.7]/sqrt(25);

% intermediate parameters
M1 = U1*V1';
K1 = V1*S1*V1';
Miu21 = mat2vec(M1);
Sigma2Inv1 = [K1,zeros(3),zeros(3)
             zeros(3),K1,zeros(3)
             zeros(3),zeros(3),K1];

% tangent space
Omega1 = skew(V1*[1;0;0]);
Omega2 = skew(V1*[0;1;0]);
Omega3 = skew(V1*[0;0;1]);
t1 = mat2vec(M1*Omega1)/sqrt(2);
t2 = mat2vec(M1*Omega2)/sqrt(2);
t3 = mat2vec(M1*Omega3)/sqrt(2);
n = null([t1';t2';t3']);
Rt = [t1';t2';t3';n'];

P1 = [PTilde1,zeros(1,6)]*Rt;

% other intermediate parameters
F1 = U1*S1*V1';
Miuc1 = @(R)Miu1+P1*Sigma2Inv1*(mat2vec(R)-Miu21);
Sigmac1 = Sigma1-P1*Sigma2Inv1*P1';

%% distribution 2
% parameters
T = expRM([0.3,0,0]);
Miu2 = Miu1;
Sigma2 = Sigma1+PTilde1*(T*0.5*diag([0,S1(1,1)-S1(2,2),S1(1,1)+S1(2,2)])*T'-...
    0.5*diag([0,S1(1,1)-S1(2,2),S1(1,1)+S1(2,2)]))*PTilde1';
U2 = U1*T;
V2 = V1*T';
S2 = S1;
PTilde2 = PTilde1*T;

% intermediate parameters
M2 = U2*V2';
K2 = V2*S2*V2';
Miu22 = mat2vec(M2);
Sigma2Inv2 = [K2,zeros(3),zeros(3)
             zeros(3),K2,zeros(3)
             zeros(3),zeros(3),K2];

% tangent space
Omega1 = skew(V2*[1;0;0]);
Omega2 = skew(V2*[0;1;0]);
Omega3 = skew(V2*[0;0;1]);
t1 = mat2vec(M2*Omega1)/sqrt(2);
t2 = mat2vec(M2*Omega2)/sqrt(2);
t3 = mat2vec(M2*Omega3)/sqrt(2);
n = null([t1';t2';t3']);
Rt = [t1';t2';t3';n'];

P2 = [PTilde2,zeros(1,6)]*Rt;

% other intermediate parameters
F2 = U2*S2*V2';
Miuc2 = @(R)Miu2+P2*Sigma2Inv2*(mat2vec(R)-Miu22);
Sigmac2 = Sigma2-P2*Sigma2Inv2*P2';

% test equivalence
FError = F1-F2;
SigmacError = Sigmac1-Sigmac2;

Nt = 100;
theta1 = linspace(-pi,pi,Nt);
theta2 = linspace(0,pi,Nt);
theta3 = linspace(-pi,pi,Nt);
MiucError = zeros(Nt^3,1);
for nt1 = 1:Nt
    for nt2 = 1:Nt
        for nt3 = 1:Nt
            index = Nt^2*(nt1-1)+Nt*(nt2-1)+nt3;
            RTest = eul2rm([theta1(nt1),theta2(nt2),theta3(nt3)]);
            MiucError(index) = Miuc1(RTest)-Miuc2(RTest);
        end
    end
end
maxMiucError = max(abs(MiucError));

end


function [ vec ] = mat2vec( mat )

vec = [mat(1,:)'; mat(2,:)'; mat(3,:)'];

end

