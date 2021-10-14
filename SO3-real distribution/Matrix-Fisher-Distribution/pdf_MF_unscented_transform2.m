function [ R, w ] = pdf_MF_unscented_transform2( F, lambda, d )

if nargin < 2 || isempty(lambda) || lambda<0 || lambda>1
    lambda = 0.1;
end

[U,S,V] = psvd(F);

if nargin < 3 || isempty(d)
    s = diag(S);
    d = pdf_MF_moment(s,false);
end

R = zeros(3,3,7);
w = zeros(7,1);

% identity
sumd = sum(d);
w(1) = lambda/3*sumd;
R(:,:,1) = U*V';

% along e1
w(2:3) = 1/6*(1+d(1)-d(2)-d(3)) + (1-lambda)/18*sumd;
theta = acos(1-(1+d(1)-d(2)-d(3))/(4*w(2)));
Q = [1,0,0;0,cos(theta),-sin(theta);0,sin(theta),cos(theta)];
R(:,:,2) = U*Q*V';
R(:,:,3) = U*Q'*V';

% along e2
w(4:5) = 1/6*(1+d(2)-d(1)-d(3)) + (1-lambda)/18*sumd;
theta = acos(1-(1+d(2)-d(1)-d(3))/(4*w(4)));
Q = [cos(theta),0,sin(theta);0,1,0;-sin(theta),0,cos(theta)];
R(:,:,4) = U*Q*V';
R(:,:,5) = U*Q'*V';

% along e1
w(6:7) = 1/6*(1+d(3)-d(1)-d(2)) + (1-lambda)/18*sumd;
theta = acos(1-(1+d(3)-d(1)-d(2))/(4*w(6)));
Q = [cos(theta),-sin(theta),0;sin(theta),cos(theta),0;0,0,1];
R(:,:,6) = U*Q*V';
R(:,:,7) = U*Q'*V';

end

