function [ V ] = wrapRotVec( V )
% wrap rotation vector to have a norm between [0,pi]
% Last edited by Weixin Wang, December 18, 2018

u = V./norm(V);
while norm(V) > pi
    V = V-2*pi*u;
end

end

