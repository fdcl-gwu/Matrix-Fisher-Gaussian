function [ pError ] = testEquivalence(  )

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
P1 = [0.5,-0.4,0.7]/sqrt(25)/3;

% density
F1 = U1*S1*V1';
Miuc1 = @(R) Miu1+P1*vee(mulRot(mulRot(U1',R,0),V1*S1,0)...
    -mulRot(mulRot(S1*V1',permute(R,[2,1,3]),0),U1,0),0,0);
Sigmac1 = Sigma1-P1*(trace(S1)*eye(3)-S1)*P1';
cF1 = pdf_MF_normal(diag(S1));
cG1 = 1/sqrt(2*pi)/sqrt(det(Sigmac1));
p1 = @(x,R) 1/cF1*exp(trace3(mulRot(repmat(F1',1,1,size(R,3)),R,0)))...
    .*(1/cG1*exp(-1/2*quadratic(Sigmac1^-1,x-Miuc1(R)))');

%% distribution 2
% parameters
T = expRot([rand,0,0]);
D3 = diag([1,1,-1]);

Miu2 = Miu1;
Sigma2 = Sigma1 + P1*(T*(trace(S1)*eye(3)-S1)*T'-(trace(S1)*eye(3)-S1))*P1';
U2 = U1*T;
V2 = V1*D3*T*D3;
S2 = S1;
P2 = P1*T;

% density
F2 = U2*S2*V2';
Miuc2 = @(R) Miu2+P2*vee(mulRot(mulRot(U2',R,0),V2*S2,0)...
    -mulRot(mulRot(S2*V2',permute(R,[2,1,3]),0),U2,0),0,0);
Sigmac2 = Sigma2-P2*(trace(S2)*eye(3)-S2)*P2';
cF2 = pdf_MF_normal(diag(S2));
cG2 = 1/sqrt(2*pi)/sqrt(det(Sigmac2));
p2 = @(x,R) 1/cF2*exp(trace3(mulRot(repmat(F2',1,1,size(R,3)),R,0)))...
    .*(1/cG2*exp(-1/2*quadratic(Sigmac2^-1,x-Miuc2(R)))');

%% test equivalence

Nt = 100;
theta1 = linspace(-pi,pi,Nt);
theta1 = repmat(theta1,Nt^2,1);
theta1 = reshape(theta1,1,[]);

theta2 = linspace(-pi,pi,Nt);
theta2 = repmat(theta2,Nt,1);
theta2 = reshape(theta2,1,[]);
theta2 = repmat(theta2,1,Nt);

theta3 = linspace(-pi,pi,Nt);
theta3 = repmat(theta3,1,Nt^2);

theta = [theta1;theta2;theta3];
R = expRot(theta);

x = linspace(-3,3,Nt);

diffp = zeros(Nt,1);
parfor nt = 1:Nt
    p1nt = p1(x(nt),R);
    p2nt = p2(x(nt),R);
    
    [diffp(nt),ind] = max(abs(p1nt-p2nt));
    if max(p1nt(ind),p2nt(ind))>0
        diffp(nt) = diffp(nt)/max(p1nt(ind),p2nt(ind));
    end
end

pError = max(diffp);

if ~any(strcmp(pathCell,getAbsPath('Matrix-Fisher-Distribution')))
    addpath('Matrix-Fisher-Distribution');
end
if ~any(strcmp(pathCell,getAbsPath('..\rotation3d')))
    addpath('..\rotation3d');
end

end


function [ result ] = quadratic( X, a )
% X is a 3-by-3 matrix, a is a 3-by-n matrix
% returns the quadratic form a'*X*a for each column of a

if size(a,1) == 1
    result = X*a.*a;
else
    N = size(a,2);
    result = zeros(N,1);
    for n = 1:N
        result(n) = a(:,n)'*X*a(:,n);
    end
end

end


function [ tr ] = trace3( F )

tr = reshape(F(1,1,:)+F(2,2,:)+F(3,3,:),[],1);

end

