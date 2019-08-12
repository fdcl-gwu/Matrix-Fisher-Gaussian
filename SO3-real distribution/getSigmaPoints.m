function [ x, R, w ] = getSigmaPoints( Miu, Sigma, PTilde, U, S, V, wM, w0 )

pathCell = regexp(path, pathsep, 'split');
if ~any(strcmp(pathCell,getAbsPath('Matrix-Fisher-Distribution')))
    addpath('Matrix-Fisher-Distribution');
end
if ~any(strcmp(pathCell,getAbsPath('..\rotation3d')))
    addpath('..\rotation3d');
end

% default values
n = size(Miu,1);
if ~exist('w0','var')
    w0 = 1/((n+3)*2+1);
end
if ~exist('wM','var')
    wM = (1-w0)/(n+3)*3;
end

% declare the size
x = zeros(n,(n+3)*2+1);
R = zeros(3,3,(n+3)*2+1);
w = zeros(1,(n+3)*2+1);

%% Matrix Fisher part
% calculate sigma
s = diag(S);
c = pdf_MF_normal(s);
dc = pdf_MF_normal_deriv(s);

cyc = {[1,2,3],[2,3,1],[3,1,2]};
wm = cell(3,1);
cost = cell(3,1);
for i = 1:3
    if s(cyc{i}(2))+s(cyc{i}(3))>=1
        cost{i} = @(sigma)sigma+(1-sigma)*(log(c)-s(cyc{i}(1)))/...
            (s(cyc{i}(2))+s(cyc{i}(3)));
    else
        cost{i} = @(sigma)(sigma+(1-sigma)*(log(c)-s(cyc{i}(1)))+1/2)*...
            (s(cyc{i}(2))+s(cyc{i}(3)))-1/2;
    end
    wm{i} = @(sigma)1/(4*(1-cost{i}(sigma)))*...
        (1/c*(dc(cyc{i}(1))-dc(cyc{i}(2))-dc(cyc{i}(3)))+1);
end

options = optimoptions('fsolve','Display','off');
sumwm = @(sigma)2*(wm{1}(sigma)+wm{2}(sigma)+wm{3}(sigma));
sigma = fsolve(@(sigma)sumwm(sigma)-wM,0,options);

% calculate sigma points and weights
fR = @(R)[trace(S*V'*R'*U*skew([1,0,0]));
    trace(S*V'*R'*U*skew([0,1,0]));
    trace(S*V'*R'*U*skew([0,0,1]))]/sqrt(2);
for i = 1:3
    theta = acos(cost{i}(sigma));
    e = zeros(3,1);
    e(i) = 1;
    R(:,:,2*n+2*i-1) = U*expRM(theta*e)*V';
    R(:,:,2*n+2*i) = U*expRM(-theta*e)*V';
    x(:,2*n+2*i-1) = Miu+PTilde*fR(R(:,:,2*n+2*i-1));
    x(:,2*n+2*i) = Miu+PTilde*fR(R(:,:,2*n+2*i));
    w([2*n+2*i-1,2*n+2*i]) = wm{i}(sigma);
end

% calculate wM again to eliminate numerical errors
wM = sum(w(2*n+1:2*n+6));

%% Gaussian part
wG = 1-w0-wM;
SigmaTildeInv = diag([s(2)+s(3),s(1)+s(3),s(1)+s(2)])/2;
Sigmac = Sigma-PTilde*SigmaTildeInv*PTilde';
for i = 1:n
    e = zeros(n,1);
    e(i) = 1;
    y = sqrt(n/wG)*e;
    x(:,2*i-1) = sqrtm(Sigmac)*y+Miu;
    x(:,2*i) = -sqrtm(Sigmac)*y+Miu;
    R(:,:,2*i-1) = U*V';
    R(:,:,2*i) = U*V';
    w([2*i-1,2*i]) = wG/2;
end

%% mean point
x(:,end) = Miu;
R(:,:,end) = U*V';
w(end) = w0;

%%
if ~any(strcmp(pathCell,getAbsPath('Matrix-Fisher-Distribution')))
    rmpath('Matrix-Fisher-Distribution');
end
if ~any(strcmp(pathCell,getAbsPath('..\rotation3d')))
    rmpath('..\rotation3d');
end

end

