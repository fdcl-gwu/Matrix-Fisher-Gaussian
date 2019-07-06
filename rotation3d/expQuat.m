function [ Q ] = expQuat( V )
% Converting rotation vector to quaternion (exponential map)
% Reference: Joan Sola, Quaternion kinematics for the error-state Kalman filter, 2017
% Last edited by Weixin Wang, December 18, 2018

if length(V) ~= 3
    error('input vector size is not 3');
end

fi = sqrt(V(1)^2+V(2)^2+V(3)^2);
if fi == 0
    Q = [1,0,0,0];
else
    u = V/fi;
    Q(1) = cos(fi/2);
    Q(2:4) = u*sin(fi/2);
end

Q = Q';

end

