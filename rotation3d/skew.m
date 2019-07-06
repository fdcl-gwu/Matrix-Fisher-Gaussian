function [ M ] = skew( V )
% skew operator for a 3-vector
% Reference: Joan Sola, Quaternion kinematics for the error-state Kalman filter, 2017
% Last edited by Weixin Wang, December 18, 2018

if length(V) ~= 3
    error('input vector size is not 3');
end

M = zeros(3);
M(1,2) = -V(3);
M(1,3) = V(2);
M(2,1) = V(3);
M(2,3) = -V(1);
M(3,1) = -V(2);
M(3,2) = V(1);

end

