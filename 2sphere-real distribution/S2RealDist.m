function [  ] = S2RealDist(  )

% parameters
n1 = 1;
n2 = 2;

Miu1 = 0;
Miu2 = [1;1;1]/sqrt(3);
Sigma = 1;
k = 10;
P = [0.5,0.5,0.5]*Sigma/sqrt(k);

% intermediate parameters
Miuc = @(x2)Miu1+k*P*(x2-Miu2);
Sigmac = Sigma-k*(P*P');

% Normalizing constant
c1 = 1/sqrt((2*pi)^n1*det(Sigmac));

p = n2+1;
c2 = (k/2)^(p/2-1)/(gamma(p/2)*besseli(p/2-1,k));

% density
f = @(x1,x2)1/c1*exp(-1/2*(x1-Miuc(x2))'*Sigmac^-1*(x1-Miuc(x2)))*...
    1/c2*exp(k*x2'*Miu2);

% grid
Nt1 = 100;
Nt2 = 100;
theta1 = linspace(0,pi,Nt1);
theta2 = linspace(-pi,pi,Nt2);
x2_1 = cos(theta2).*sin(theta1)';
x2_2 = sin(theta2).*sin(theta1)';
x2_3 = repmat(cos(theta1)',1,Nt1);

Nx1 = 11;
x1 = linspace(-1,1,Nx1);

% colormap
color = zeros(Nt1,Nt2,Nx1);
for nt1 = 1:Nt1
    for nt2 = 1:Nt2
        for nx1 = 1:Nx1
            color(nt1,nt2,nx1) = f(x1(nx1),[x2_1(nt1,nt2);x2_2(nt1,nt2);x2_3(nt1,nt2)]);
        end
    end
end
cmax = max(max(max(color)));

% plot
for nx1 = 1:Nx1
    figure;
    surf(x2_1,x2_2,x2_3,color(:,:,nx1),'LineStyle','none');
    caxis([0,cmax]);
    axis equal;
    view([1,1,0.2]);
end

end

