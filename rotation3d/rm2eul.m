function [ eul ] = rm2eul( rm )
% converting rotation matrix to Euler angles (3-2-1 body)
% Reference: en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
% Last edited by Weixin Wang, December 18, 2018

n = size(rm,3);
if n > 1
    rm11 = reshape(rm(1,1,:),n,1);
    rm21 = reshape(rm(2,1,:),n,1);
    rm31 = reshape(rm(3,1,:),n,1);
    rm32 = reshape(rm(3,2,:),n,1);
    rm33 = reshape(rm(3,3,:),n,1);
else
    rm11 = rm(1,1);
    rm21 = rm(2,1);
    rm31 = rm(3,1);
    rm32 = rm(3,2);
    rm33 = rm(3,3);
end

eul(:,2) = -asin(rm31);
eul(:,1) = atan2(rm32./cos(eul(:,2)),rm33./cos(eul(:,2)));
eul(:,3) = atan2(rm21./cos(eul(:,2)),rm11./cos(eul(:,2)));

end

