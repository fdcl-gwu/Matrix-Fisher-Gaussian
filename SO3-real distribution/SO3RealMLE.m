function [ Miu, Sigma, PTilde, UO, SO, VO, invariant1, invariant2 ] = SO3RealMLE( x, R )
% let x be N-by-Ns, R be 3-by-3-by-Ns

N = size(x,1);
Ns = size(R,3);

% empirical moments
Ex = sum(x,2)/Ns;
ER = sum(R,3)/Ns;

covxx = zeros(N,N);
for ns = 1:Ns
    covxx = covxx+(x(:,ns)-Ex)*(x(:,ns)-Ex)'/Ns;
end

% Matrix Fisher part
[UO,DO,VO] = usvd(ER,true);
SO = diag(pdf_MF_M2S(diag(DO)));

% tangent space
MO = UO*VO';
Omega1 = skew(VO*[1;0;0]);
Omega2 = skew(VO*[0;1;0]);
Omega3 = skew(VO*[0;0;1]);
t1 = mat2vec(MO*Omega1)/sqrt(2);
t2 = mat2vec(MO*Omega2)/sqrt(2);
t3 = mat2vec(MO*Omega3)/sqrt(2);
n = null([t1';t2';t3']);
Rt = [t1';t2';t3';n'];

% f(eta)
fEta = zeros(3,Ns);
for ns = 1:Ns
    expEta = UO'*R(:,:,ns)*VO;
    fEta(1,ns) = trace(SO*expEta'*skew([1,0,0]))/sqrt(2);
    fEta(2,ns) = trace(SO*expEta'*skew([0,1,0]))/sqrt(2);
    fEta(3,ns) = trace(SO*expEta'*skew([0,0,1]))/sqrt(2);
end

% other empirical moments
EfEta = sum(fEta,2)/Ns;

covxfEta = zeros(N,N);
covfEtafEta = zeros(N,N);
for ns = 1:Ns
    covxfEta = covxfEta+(x(:,ns)-Ex)*(fEta(:,ns)-EfEta)'/Ns;
    covfEtafEta = covfEtafEta+(fEta(:,ns)-EfEta)*(fEta(:,ns)-EfEta)'/Ns;
end

% correlation part
PTilde = covxfEta*covfEtafEta^-1;
P = [PTilde,zeros(1,6)]*Rt;

% Gaussian part
SigmaTilde2Inv = diag([SO(2,2)+SO(3,3),SO(1,1)+SO(3,3),SO(1,1)+SO(2,2)])/2;
Miu = Ex-PTilde*EfEta;
Sigma = covxx-covxfEta*covfEtafEta^-1*covxfEta'+...
    PTilde*SigmaTilde2Inv*PTilde';

% invariants
K = VO*SO*VO';
Sigma2Inv = [K,zeros(3),zeros(3)
    zeros(3),K,zeros(3)
    zeros(3),zeros(3),K];
invariant1 = Sigma-P*Sigma2Inv*P';
invariant2 = P*Sigma2Inv;

end


function [ vec ] = mat2vec( mat )

vec = [mat(1,:)'; mat(2,:)'; mat(3,:)'];

end

