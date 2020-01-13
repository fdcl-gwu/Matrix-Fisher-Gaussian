function [ R, w ] = MFGetSigmaPoints( F, w0 )

filePath = mfilename('fullpath');
pathCell = regexp(path, pathsep, 'split');
if ~any(strcmp(pathCell,getAbsPath('..\Matrix-Fisher-Distribution',filePath)))
    addpath(getAbsPath('..\Matrix-Fisher-Distribution',filePath));
end
if ~any(strcmp(pathCell,getAbsPath('..\..\rotation3d',filePath)))
    addpath(getAbsPath('..\..\rotation3d',filePath));
end

% default values
if ~exist('w0','var') || isempty(w0)
    w0 = 1/7;
end

% calculate sigma (use scaled normalization)
[U,S,V] = psvd(F);
s = diag(S);
c = pdf_MF_normal(s,true);
dc = pdf_MF_normal_deriv(s,false,true);

cyc = {[1,2,3],[2,3,1],[3,1,2]};
wm = cell(3,1);
cost = cell(3,1);
for i = 1:3
    if s(cyc{i}(2))+s(cyc{i}(3))>=1
        cost{i} = @(sigma)sigma+(1-sigma)*(log(c)+sum(s)-s(cyc{i}(1)))/...
            (s(cyc{i}(2))+s(cyc{i}(3)));
    else
        cost{i} = @(sigma)(sigma+(1-sigma)*(log(c)+sum(s)-s(cyc{i}(1)))+1/2)*...
            (s(cyc{i}(2))+s(cyc{i}(3)))-1/2;
    end
    wm{i} = @(sigma)1/(4*(1-cost{i}(sigma)))*...
        (1/c*(dc(cyc{i}(1))-dc(cyc{i}(2))-dc(cyc{i}(3))));
end

% minimum of sigma {cos(theta_i)>-1/2}
options = optimoptions('fsolve','Display','off');
sigmaMin = max(cellfun(@(eqn)fsolve(@(sigma)eqn(sigma)+1/2,0,options),cost));

% calculate sigma
sumwm = @(sigma)2*(wm{1}(sigma)+wm{2}(sigma)+wm{3}(sigma));
sigma = fsolve(@(sigma)sumwm(sigma)-(1-w0),0,options);
if sigma < sigmaMin
    sigma = sigmaMin;
end

% calculate sigma points and weights
w = zeros(1,7);
R = zeros(3,3,7);
for i = 1:3
    theta = acos(cost{i}(sigma));
    e = zeros(3,1);
    e(i) = 1;
    R(:,:,2*i-1) = U*expRot(theta*e)*V';
    R(:,:,2*i) = U*expRot(-theta*e)*V';
    w([2*i-1,2*i]) = wm{i}(sigma);
end

% the point at mean
w(end) = 1-sum(w(1:6));
R(:,:,end) = U*V';

if ~any(strcmp(pathCell,getAbsPath('..\Matrix-Fisher-Distribution',filePath)))
    rmpath(getAbsPath('..\Matrix-Fisher-Distribution',filePath));
end
if ~any(strcmp(pathCell,getAbsPath('..\..\rotation3d',filePath)))
    rmpath(getAbsPath('..\..\rotation3d',filePath));
end

end

