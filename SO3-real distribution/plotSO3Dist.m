function [ c ] = plotSO3Dist( theta1, theta2, f )
% f is the density function handle
% s1, s2, s3 are coordinates on the ball

addpath('..\rotation3d');

Nt1 = length(theta1);
Nt2 = length(theta2);

% grid of rotation about the axis of interest
Nt3 = 50;
theta3 = linspace(-pi,pi,Nt3);

% marginal density for each column the of rotation matrix
f1 = zeros(Nt1,Nt2);
f2 = zeros(Nt1,Nt2);
f3 = zeros(Nt1,Nt2);
tic;
parfor nt1 = 1:Nt1
    for nt2 = 1:Nt2
        axis = [cos(theta1(nt1))*sin(theta2(nt2));...
            sin(theta1(nt1))*sin(theta2(nt2)); cos(theta2(nt2))];
        RAxis = zeros(3,3,Nt3);
        for nt3 = 1:Nt3
            RAxis(:,:,nt3) = expRM(axis*theta3(nt3));
        end
        
        ortAxes = null(axis');
        ortAxisGrid = zeros(3,2,Nt3);
        for nt3 = 1:Nt3
            ortAxisGrid(:,:,nt3) = RAxis(:,:,nt3)*ortAxes;
        end
        
        R1 = zeros(3,3,Nt3);
        R2 = zeros(3,3,Nt3);
        R3 = zeros(3,3,Nt3);
        for nt3 = 1:Nt3
            R1(:,:,nt3) = [axis,ortAxisGrid(:,:,nt3)];
            if det(R1(:,:,nt3))<0
                R1(:,:,nt3) = [axis,flip(ortAxisGrid(:,:,nt3),2)];
            end
            
            R2(:,:,nt3) = [ortAxisGrid(:,1,nt3),axis,...
                ortAxisGrid(:,2,nt3)];
            if det(R2(:,:,nt3))<0
                R2(:,:,nt3) = [ortAxisGrid(:,2,nt3),axis,...
                    ortAxisGrid(:,1,nt3)];
            end
            
            R3(:,:,nt3) = [ortAxisGrid(:,:,nt3),axis];
            if det(R3(:,:,nt3))<0
                R3(:,:,nt3) = [flip(ortAxisGrid(:,:,nt3),2),axis];
            end
        end
        
        for nt3 = 1:Nt3
            f1(nt1,nt2) = f1(nt1,nt2)+f(R1(:,:,nt3))*(2*pi)/Nt3;
            f2(nt1,nt2) = f2(nt1,nt2)+f(R2(:,:,nt3))*(2*pi)/Nt3;
            f3(nt1,nt2) = f3(nt1,nt2)+f(R3(:,:,nt3))*(2*pi)/Nt3;
        end
    end
end
toc;

% colormap
c = f1+f2+f3;

rmpath('..\rotation3d');

end

