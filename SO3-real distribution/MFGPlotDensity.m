function [  ] = MFGPlotDensity(  )

pathCell = regexp(path, pathsep, 'split');
if ~any(strcmp(pathCell,getAbsPath('Matrix-Fisher-Distribution')))
    addpath('Matrix-Fisher-Distribution');
end
if ~any(strcmp(pathCell,getAbsPath('..\rotation3d')))
    addpath('..\rotation3d');
end

% parameters
n1 = 1;

Miu = 0;
Sigma = 1;
U = expRot([0,0,0]);
V = expRot([0,0,0]);
S = diag([25,25,25]);
F = U*S*V';
P = [0,0,0.7]/sqrt(25+25);

% intermediate parameters
[U,S,V] = psvd(F);
Sigma2Inv = trace(S)*eye(3)-S;

% other intermediate parameters
Q = @(R)U'*R*V;
Miuc = @(R)Miu+P*vee(Q(R)*S-S*Q(R)');
Sigmac = Sigma-P*Sigma2Inv*P';

% Normalizing constant
c1 = 1/sqrt((2*pi)^n1*det(Sigmac));
c2 = pdf_MF_normal(diag(S));

% density
f = @(x1,R)1/c1*exp(-1/2*(x1-Miuc(R))'*Sigmac^-1*(x1-Miuc(R)))*...
    1/c2*exp(trace(F*R'));

% spherical grid
Nt1 = 100;
Nt2 = 100;
theta1 = linspace(-pi,pi,Nt2);
theta2 = linspace(0,pi,Nt1);
s1 = cos(theta1)'.*sin(theta2);
s2 = sin(theta1)'.*sin(theta2);
s3 = repmat(cos(theta2),Nt1,1);

% linear grid
Nx1 = 11;
x1 = linspace(-1,1,Nx1);

% color map
c = zeros(Nt1,Nt2,Nx1);
for nx1 = 1:Nx1
    fRot = @(R)f(x1(nx1),R);
    c(:,:,nx1) = plotSO3Dist(theta1,theta2,fRot);
end
cmax = max(max(max(c)));

% plot Matrix Fisher
M = U*V';
for nx1 = 1:Nx1
    figure; hold on;
    surf(s1,s2,s3,c(:,:,nx1),'LineStyle','none');
    plot3([0,M(1,1)*1.1],[0,M(2,1)*1.1],[0,M(3,1)*1.1],'b');
    plot3([0,M(1,2)*1.1],[0,M(2,2)*1.1],[0,M(3,2)*1.1],'r');
    plot3([0,M(1,3)*1.1],[0,M(2,3)*1.1],[0,M(3,3)*1.1],'y');
    
    axis equal;
    view([1,1,1]);
    caxis([0,cmax]);
end

% linear grid
Nx1 = 1001;
x1 = linspace(-4,4,Nx1);

% rotation grid
Nt1 = 11;
theta1 = linspace(-pi/8,pi/8,Nt1);

% plot Gaussian
fLinear = zeros(Nx1,Nt1,3);
MiuLinear = zeros(Nt1,3);
for i = 1:3
    n = zeros(3,1);
    for nt1 = 1:Nt1
        n(i) = theta1(nt1);
        MiuLinear(nt1,i) = Miuc(M*expRot(V*n));
        for nx1 = 1:Nx1
            fLinear(nx1,nt1,i) = f(x1(nx1),M*expRot(V*n));
        end
    end
end
fmax = max(max(max(fLinear)));

for i = 1:3
    for nt1 = 1:Nt1
        figure; hold on;
        plot(x1,fLinear(:,nt1,i));
        
        xlim([-4,4]);
        ylim([0,fmax]);
        
        str = strcat('$\mu=$',num2str(MiuLinear(nt1,i)),', $\sigma=$',num2str(Sigmac));
        annotation('textbox','String',str,'Interpreter','latex');
    end
end


rmpath('Matrix-Fisher-Distribution');
rmpath('..\rotation3d');

end

