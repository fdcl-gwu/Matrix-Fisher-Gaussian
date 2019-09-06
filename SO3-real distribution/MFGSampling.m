function [ x, R ] = MFGSampling( Miu, Sigma, P, U, S, V, Ns )

pathCell = regexp(path, pathsep, 'split');
if ~any(strcmp(pathCell,getAbsPath('Matrix-Fisher-Distribution')))
    addpath('Matrix-Fisher-Distribution');
end
if ~any(strcmp(pathCell,getAbsPath('..\rotation3d')))
    addpath('..\rotation3d');
end

N = length(Miu);

% sample from Matrix Fisher distribution
F = U*S*V';
R = pdf_MF_sampling(F,Ns);

% sample from Gaussian ditribution
Sigma2Inv = diag([S(2,2)+S(3,3),S(1,1)+S(3,3),S(1,1)+S(2,2)]);
Sigmac = Sigma-P*Sigma2Inv*P';
x = zeros(N,Ns);
for ns = 1:Ns
    Q = U.'*R(:,:,ns)*V;
    fR = vee(Q*S-S*Q.');
    Miuc = Miu+P*fR;
    x(:,ns) = mvnrnd(Miuc,Sigmac);
end

end

