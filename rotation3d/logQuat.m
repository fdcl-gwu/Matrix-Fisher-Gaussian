function [ V ] = logQuat( Q )
% Converting quaternion to rotation vector (logarithm map)
% Reference: Joan Sola, Quaternion kinematics for the error-state Kalman filter, 2017
% Last edited by Weixin Wang, December 18, 2018

if Q(1) == 1
    V = [0,0,0];
    return;
end

phi = 2*atan2(norm(Q(2:4)),Q(1));
u = Q(2:4)/norm(Q(2:4));

V = (phi*u)';

end

