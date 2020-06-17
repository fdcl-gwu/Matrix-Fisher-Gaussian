function [ R, w ] = MFGetSigmaPoints( F, w0 )

if ~exist('w0','var') || isempty(w0)
    w0 = 1/7;
end

% calculate sigma (use scaled normalization)
[U,S,V] = psvd(F);
s = diag(S);
[c,dc] = pdf_MF_normal(s,1,1);

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

% minimum of sigma {cos(theta_i)>-sqrt(3)/2}
options = optimoptions('fsolve','Display','off');
sigmaMin = max(cellfun(@(eqn)fsolve(@(sigma)eqn(sigma)+sqrt(3)/2,0,options),cost));

% calculate sigma
sumwm = @(sigma)2*(wm{1}(sigma)+wm{2}(sigma)+wm{3}(sigma));
sigma = fsolve(@(sigma)sumwm(sigma)-(1-w0),0,options);
if sigma <= sigmaMin
    sigma = sigmaMin;
    wM = sumwm(sigma);
    w0 = 1-wM;
end

% calculate sigma points and weights
R = zeros(3,3,7);
w = zeros(1,7);

for i = 1:3
    theta = acos(cost{i}(sigma));
    e = zeros(3,1);
    e(i) = 1;
    R(:,:,2*i) = U*expRot(theta*e)*V';
    R(:,:,2*i+1) = U*expRot(-theta*e)*V';
    w([2*i,2*i+1]) = wm{i}(sigma);
end

R(:,:,7) = eye(3);
w(7) = w0;

end

