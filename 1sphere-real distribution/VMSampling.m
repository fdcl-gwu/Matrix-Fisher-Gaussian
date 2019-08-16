function [ theta ] = VMSampling( k, miu, Ns )

U = rand(Ns,3);

a = 1+(1+4*k^2)^(1/2);
b = (a-(2*a)^(1/2))/(2*k);
r = (1+b^2)/(2*b);

theta = zeros(1,Ns);
for n = 1:Ns
    z = cos(pi*U(n,1));
    f = (1+r*z)/(r+z);
    c = k*(r-f);
    
    while log(c/U(n,2))+1-c < 0
        if c*(2-c)-U(n,2)>0
            break;
        else
            U(n,:) = rand(1,3);
            z = cos(pi*U(n,1));
            f = (1+r*z)/(r+z);
            c = k*(r-f);
        end
    end
    
    theta(n) = miu+sign(U(n,3)-0.5)*acos(f);
end

end

