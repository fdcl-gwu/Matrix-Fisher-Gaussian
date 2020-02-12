function [ x, R ] = MFGSampling( Miu, Sigma, P, U, S, V, Ns )

pathCell = regexp(path, pathsep, 'split');
if ~any(strcmp(pathCell,getAbsPath('Matrix-Fisher-Distribution')))
    addpath('Matrix-Fisher-Distribution');
end
if ~any(strcmp(pathCell,getAbsPath('..\rotation3d')))
    addpath('..\rotation3d');
end

N = size(Miu,1);

% sample from canonical MFG
y = mvnrnd(zeros(N,1),eye(N),Ns)';
Q = pdf_MF_sampling(S,Ns);

% transfomr back to MFG
SigmaMInv = diag([S(2,2)+S(3,3),S(1,1)+S(3,3),S(1,1)+S(2,2)]);
Sigmac = Sigma-P*SigmaMInv*P';

R = zeros(3,3,Ns);
x = zeros(3,Ns);
for ns = 1:Ns
    R(:,:,ns) = U*Q(:,:,ns)*V';
    gR = vee(Q(:,:,ns)*S-S*Q(:,:,ns)');
    x(:,ns) = sqrtm(Sigmac)*y(:,ns)+Miu+P*gR;
end

if ~any(strcmp(pathCell,getAbsPath('Matrix-Fisher-Distribution')))
    rmpath('Matrix-Fisher-Distribution');
end
if ~any(strcmp(pathCell,getAbsPath('..\rotation3d')))
    rmpath('..\rotation3d');
end

end

