function [ Q ]=rot2qua( R )
% converting rotation matrix to quaternion
% Reference: www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/
% Last edited by Weixin Wang, December 18, 2018

tr = R(1,1)+R(2,2)+R(3,3);
if tr > 0
    S = sqrt(tr+1)*2;
    qr = S/4;
    qx = (R(3,2)-R(2,3))/S;
    qy = (R(1,3)-R(3,1))/S;
    qz = (R(2,1)-R(1,2))/S;
elseif R(1,1)>R(2,2) && R(1,1)>R(3,3)
    S = sqrt(1+R(1,1)-R(2,2)-R(3,3))*2;
    qr = (R(3,2)-R(2,3))/S;
    qx = S/4;
    qy = (R(1,2)+R(2,1))/S;
    qz = (R(1,3)+R(3,1))/S;
elseif R(2,2)>R(3,3)
    S = sqrt(1+R(2,2)-R(1,1)-R(3,3))*2;
    qr = (R(1,3)-R(3,1))/S;
    qx = (R(1,2)+R(2,1))/S;
    qy = S/4;
    qz = (R(2,3)+R(3,2))/S;
else
    S = sqrt(1+R(3,3)-R(1,1)-R(2,2))*2;
    qr = (R(2,1)-R(1,2))/S;
    qx = (R(1,3)+R(3,1))/S;
    qy = (R(2,3)+R(3,2))/S;
    qz = S/4;
end

Q = [qr,qx,qy,qz];
