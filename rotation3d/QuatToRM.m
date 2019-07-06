function [ RM ] = QuatToRM( Q )
% converting quaternion to rotation matrix
% Reference: Joan Sola, Quaternion kinematics for the error-state Kalman filter, 2017
% Last edited by Weixin Wang, December 18, 2018

a = Q(1);b=Q(2);c=Q(3);d=Q(4);

RM = zeros(3);
RM(1,1) = a^2+b^2-c^2-d^2;
RM(1,2) = 2*b*c-2*a*d;
RM(1,3) = 2*b*d+2*a*c;
RM(2,1) = 2*b*c+2*a*d;
RM(2,2) = a^2-b^2+c^2-d^2;
RM(2,3) = 2*c*d-2*a*b;
RM(3,1) = 2*b*d-2*a*c;
RM(3,2) = 2*c*d+2*a*b;
RM(3,3) = a^2-b^2-c^2+d^2;

end

