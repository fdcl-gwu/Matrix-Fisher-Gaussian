function [ V ] = logRM( RM )
% Converting rotation matrix to rotation vector (logarithm map)
% Reference: Joan Sola, Quaternion kinematics for the error-state Kalman filter, 2017
% Last edited by Weixin Wang, December 18, 2018

phi = acos((trace(RM)-1)/2);
u = 1/(2*sin(phi))*[RM(3,2)-RM(2,3);RM(1,3)-RM(3,1);RM(2,1)-RM(1,2)];

V = phi*u;

end

