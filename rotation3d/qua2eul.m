% Author: Lauro Ojeda, 2003-2015
% Bibliography:
% David Titterton and John Weston, Strapdown Inertial Navigation Technology, 2nd Edition
% pp 45
% Converts quaternion to Euler representation
function euler = qua2eul(quaternion)

if size(quaternion,1) == 4
    quaternion = quaternion';
end

a = quaternion(:,1);
b = quaternion(:,2);
c = quaternion(:,3);
d = quaternion(:,4);

% deal with gimbal lock numerical issue
sinthe = 2*(a.*c - d.*b);
sinthe(sinthe < -1) = -1;
sinthe(sinthe > 1) = 1;

phi = atan2((2*(a.*b + c.*d)),(a.^2 - b.^2 - c.^2 + d.^2));
the = asin(sinthe);
psi = atan2((2*(a.*d + b.*c)),(a.^2 + b.^2 - c.^2 - d.^2));
euler = [phi,the,psi];
