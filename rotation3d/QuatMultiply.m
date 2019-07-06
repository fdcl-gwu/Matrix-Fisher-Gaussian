function [ qout ] = QuatMultiply( q1, q2 )
% calculate the quaternion multiplication of q1 x q2
% Reference: Joan Sola, Quaternion kinematics for the error-state Kalman filter, 2017
% Last edited by Weixin Wang, December 18, 2018

qout = zeros(4,1);
qout(1) = q1(1)*q2(1)-q1(2)*q2(2)-q1(3)*q2(3)-q1(4)*q2(4);
qout(2) = q1(1)*q2(2)+q1(2)*q2(1)+q1(3)*q2(4)-q1(4)*q2(3);
qout(3) = q1(1)*q2(3)-q1(2)*q2(4)+q1(3)*q2(1)+q1(4)*q2(2);
qout(4) = q1(1)*q2(4)+q1(2)*q2(3)-q1(3)*q2(2)+q1(4)*q2(1);

end

