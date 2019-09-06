function [ Miu, Sigma, P, U, S, V ] = MFGMLEAppro( x, R, w )
% let x be N-by-Ns, R be 3-by-3-by-Ns

filePath = mfilename('fullpath');
pathCell = regexp(path, pathsep, 'split');
if ~any(strcmp(pathCell,getAbsPath('Matrix-Fisher-Distribution',filePath)))
    addpath(getAbsPath('Matrix-Fisher-Distribution',filePath));
end
if ~any(strcmp(pathCell,getAbsPath('..\rotation3d',filePath)))
    addpath(getAbsPath('..\rotation3d',filePath));
end

N = size(x,1);
Ns = size(R,3);

if ~exist('w','var')
    w = ones(1,Ns)/Ns;
end

% empirical moments
Ex = sum(x.*w,2);
ER = sum(R.*reshape(w,1,1,[]),3);

covxx = zeros(N,N);
for ns = 1:Ns
    covxx = covxx+(x(:,ns)-Ex)*(x(:,ns)-Ex)'*w(ns);
end

% Matrix Fisher part
[Up,Dp,Vp] = usvd(ER,true);
Sp = diag(pdf_MF_M2S(diag(Dp)));
[U,S,V] = usvd(Up*Sp*Vp');

% f(eta)
fEta = zeros(3,Ns);
for ns = 1:Ns
    expEta = U'*R(:,:,ns)*V;
    fEta(1,ns) = trace(S*expEta'*hat([1,0,0]));
    fEta(2,ns) = trace(S*expEta'*hat([0,1,0]));
    fEta(3,ns) = trace(S*expEta'*hat([0,0,1]));
end

% other empirical moments
EfEta = sum(fEta.*w,2);

covxfEta = zeros(N,N);
covfEtafEta = zeros(N,N);
for ns = 1:Ns
    covxfEta = covxfEta+(x(:,ns)-Ex)*(fEta(:,ns)-EfEta)'*w(ns);
    covfEtafEta = covfEtafEta+(fEta(:,ns)-EfEta)*(fEta(:,ns)-EfEta)'*w(ns);
end

% correlation part
P = covxfEta*covfEtafEta^-1;

% Gaussian part
SigmaTilde2Inv = diag([S(2,2)+S(3,3),S(1,1)+S(3,3),S(1,1)+S(2,2)]);
Miu = Ex-P*EfEta;
Sigma = covxx-covxfEta*covfEtafEta^-1*covxfEta'+...
    P*SigmaTilde2Inv*P';

if ~any(strcmp(pathCell,getAbsPath('Matrix-Fisher-Distribution',filePath)))
    rmpath(getAbsPath('Matrix-Fisher-Distribution',filePath));
end
if ~any(strcmp(pathCell,getAbsPath('..\rotation3d',filePath)))
    rmpath(getAbsPath('..\rotation3d',filePath));
end

end


function [ vec ] = mat2vec( mat )

vec = [mat(1,:)'; mat(2,:)'; mat(3,:)'];

end

