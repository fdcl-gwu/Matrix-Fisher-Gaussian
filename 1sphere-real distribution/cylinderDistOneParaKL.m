function [ div ] = cylinderDistOneParaKL( miu1, sigma1, P1, k1, theta01,...
    miu2, sigma2, P2, k2, theta02)

% Von Mises part
div = -log(besseli(0,k1)/besseli(0,k2));
A = @(k)besseli(1,k)/besseli(0,k);
div = div+k1*A(k1)-k2*A(k1)*cos(theta01-theta02);

% Gaussian part
A2 = @(k)besseli(2,k)/besseli(0,k);
sigmac1 = sqrt(sigma1^2-k1*P1^2);
sigmac2 = sqrt(sigma2^2-k2*P2^2);
div = div-log(sigmac1/sigmac2)-1/2;
div = div+1/2*sigmac1^2/sigmac2^2+(miu1-miu2)^2/(2*sigmac2^2)...
    -k2*P2*A(k1)/sigmac2^2*(miu1-miu2)*sin(theta01-theta02)...
    +(k1^2*P1^2-2*k1*P1*k2*P2*cos(theta01-theta02)+k2^2*P2^2)/(4*sigmac2^2)...
    -A2(k1)/(4*sigmac2^2)*((k1*P1-k2*P2*cos(theta01-theta02))^2-k2^2*P2^2*sin(theta01-theta02)^2);

end

