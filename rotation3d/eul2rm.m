function [ rm ] = eul2rm( eul )
% converting Euler angles (3-2-1 body) to rotation matrix
% Reference: en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
% Last edited by Weixin Wang, December 18, 2018

rm(1,1) = cos(eul(2))*cos(eul(3));
rm(1,2) = -cos(eul(1))*sin(eul(3)) + sin(eul(1))*sin(eul(2))*cos(eul(3));
rm(1,3) = sin(eul(1))*sin(eul(3)) + cos(eul(1))*sin(eul(2))*cos(eul(3));
rm(2,1) = cos(eul(2))*sin(eul(3));
rm(2,2) = cos(eul(1))*cos(eul(3)) + sin(eul(1))*sin(eul(2))*sin(eul(3));
rm(2,3) = -sin(eul(1))*cos(eul(3)) + cos(eul(1))*sin(eul(2))*sin(eul(3));
rm(3,1) = -sin(eul(2));
rm(3,2) = sin(eul(1))*cos(eul(2));
rm(3,3) = cos(eul(1))*cos(eul(2));

end

