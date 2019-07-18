function [ Miu, Sigma, P, M, K, PTilde ] = SO3RealMLE( x, R )
% let x be N-by-Ns, R be 3-by-3-by-Ns

N = size(x,1);
Ns = size(R,3);

% empirical moments
Ex = sum(x,2)/Ns;
ER = sum(R,3)/Ns;

covxx = zeros(N,N);
covxR = zeros(N,9);
covRR = zeros(9,9);
for ns = 1:Ns
    covxx = covxx+(x(:,ns)-Ex)*(x(:,ns)-Ex)'/Ns;
    covxR = covxR+(x(:,ns)-Ex)*(mat2vec(R(:,:,ns))-mat2vec(ER))'/Ns;
    covRR = covRR+(mat2vec(R(:,:,ns))-mat2vec(ER))...
        *(mat2vec(R(:,:,ns))-mat2vec(ER))'/Ns;
end

% Matrix Fisher part
[Us,Ds,Vs] = psvd(ER);
Ss = pdf_MF_M2S(diag(Ds));
M = Us*Vs';
K = Vs*diag(Ss)*Vs';

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

% correlation part
Sigma2Inv = [K,zeros(3),zeros(3)
             zeros(3),K,zeros(3)
             zeros(3),zeros(3),K];
P = covxR*covRR^-1*Sigma2Inv^-1;
PTilde = P*Rt';

% Gaussian part
Miu = Ex+covxR*covRR^-1*(mat2vec(M)-mat2vec(ER));
Sigma = covxx-covxR*covRR^-1*covxR'+P*Sigma2Inv*P';

end


function [ vec ] = mat2vec( mat )

vec = [mat(1,:)'; mat(2,:)'; mat(3,:)'];

end

