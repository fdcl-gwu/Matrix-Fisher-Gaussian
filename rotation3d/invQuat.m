function [ invQ ] = invQuat( Q )
% inverse of a quaternion
% Reference: Joan Sola, Quaternion kinematics for the error-state Kalman filter, 2017
% Last edited by Weixin Wang, December 18, 2018

invQ = [Q(1);-Q(2);-Q(3);-Q(4)]/norm(Q)^2;

end

