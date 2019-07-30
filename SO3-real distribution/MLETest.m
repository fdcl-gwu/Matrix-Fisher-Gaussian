function [ MiuError, maxMiucError, SigmacError ] = MLETest(  )

pathCell = regexp(path, pathsep, 'split');
if ~any(strcmp(pathCell,getAbsPath('Matrix-Fisher-Distribution')))
    addpath('Matrix-Fisher-Distribution');
end
if ~any(strcmp(pathCell,getAbsPath('..\rotation3d')))
    addpath('..\rotation3d');
end

% run several times
Nsim = 100;
MiuError = zeros(Nsim,1);
maxMiucError = zeros(Nsim,1);
SigmacError = zeros(Nsim,1);

% a grid over SO(3)
Nt = 100;
theta1 = linspace(-pi,pi,Nt);
theta2 = linspace(0,pi,Nt);
theta3 = linspace(-pi,pi,Nt);
RTest = zeros(3,3,Nt^3);
for nt1 = 1:Nt
    for nt2 = 1:Nt
        for nt3 = 1:Nt
            index = Nt^2*(nt1-1)+Nt*(nt2-1)+nt3;
            RTest(:,:,index) = eul2rm([theta1(nt1),theta2(nt2),theta3(nt3)]);
        end
    end
end

% parameters
N = 1;
Ns = 100000;

Miu = 0;
Sigma = 1;
USO = expRM([0.2,0.3,0.4]);
VSO = expRM([-0.1,0.5,0.7]);
SSO = diag([25,25,-25]);
PTilde = [0.7,0.5,-0.4]/sqrt(25);

% intermediate parameters
F = USO*SSO*VSO';
[UO,SO,VO] = usvd(F,true);
MO = UO*VO';
KO = VO*SO*VO';
Miu2 = mat2vec(MO);
Sigma2Inv = [KO,zeros(3),zeros(3)
             zeros(3),KO,zeros(3)
             zeros(3),zeros(3),KO];

% tangent space
Omega1 = skew(VO*[1;0;0]);
Omega2 = skew(VO*[0;1;0]);
Omega3 = skew(VO*[0;0;1]);
t1 = mat2vec(MO*Omega1)/sqrt(2);
t2 = mat2vec(MO*Omega2)/sqrt(2);
t3 = mat2vec(MO*Omega3)/sqrt(2);
n = null([t1';t2';t3']);
Rt = [t1';t2';t3';n'];

P = [PTilde,zeros(1,6)]*Rt;

for nsim = 1:Nsim
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
    [MiuHat,SigmaHat,PTildeHat,UHat,SHat,VHat,in1Hat,in2Hat] = SO3RealMLE(xs,Rs);
    SigmacError(nsim) = in1Hat-Sigmac;
    MiuError(nsim) = Miu-MiuHat;
    MiucError = zeros(Nt^3,1);
    for nt = 1:Nt^3
        MiucError(nt) = in2Hat*mat2vec(RTest(:,:,nt))-...
            P*Sigma2Inv*mat2vec(RTest(:,:,nt));
    end
    maxMiucError(nsim) = max(abs(MiucError));
end

rmpath('Matrix-Fisher-Distribution');
rmpath('..\rotation3d');

end


function [ vec ] = mat2vec( mat )

vec = [mat(1,:)'; mat(2,:)'; mat(3,:)'];

end

