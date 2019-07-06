function [ R ] = expRM( V )
% Converting rotation vector to rotation matrix (exponential map)
% Reference: Joan Sola, Quaternion kinematics for the error-state Kalman filter, 2017
% Last edited by Weixin Wang, December 18, 2018

if length(V) ~= 3
    error('input vector size is not 3');
end

fi = norm(V);
if fi == 0
    R = eye(3);
else
    u = V/fi;
    uSkew = skew(u);
    R = eye(3) + sin(fi)*uSkew + (1-cos(fi))*uSkew^2;
end

end

