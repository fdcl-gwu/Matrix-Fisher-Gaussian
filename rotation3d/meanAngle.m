function [ meanx, meanl ] = meanAngle( x )
% Calculating the mean direction and mean resultant length of a set of angluar values
% Reference: Kanti V. Mardia, Peter E. Jupp, Directional Statistics, 2000
% Last edited by Weixin Wang, December 18, 2018

sinx = sin(x);
cosx = cos(x);
ms = mean(sinx);
mc = mean(cosx);
meanx = atan2(ms,mc);
meanl = sqrt(ms^2+mc^2);

end

