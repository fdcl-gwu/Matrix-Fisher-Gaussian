function [s FVAL NITER]=pdf_MF_M2S(d,s0)
%pdf_MF_M2S: transforms the first moments into the proper singular values
%   s=pdf_MF_M2S(d,s0) numerically solves the following equations for s
%
%       \frac{1}{c(S)}\frac{\partial c(S)}{s_i} - d_i = 0,
%
%   to find the proper singular values of the matrix Fisher distribution
%   whose first moments match with d. It is solved via Newton-Armijo
%   iteration with a polynomial line search, and the optional input variable 
%   s0 specifis the initial guess of the iteration.
%
%   [s, FVAL, niter]=pdf_MF_M2S(d,s0) returns the value of the equation at
%   s, and the number of iterations
%
%   See T. Lee, "Bayesian Attitude Estimation with the Matrix Fisher
%   Distribution on SO(3)", 2017, http://arxiv.org/abs/1710.03746
%   C. T. Kelly, "Iterative Methods For Linear And Nonlinear Equations,"
%   SIAM, 1995, Section 8.3.1.
%
%   See also PDF_MF_M2S_APPROX

bool_scaled=1;
if nargin < 2
    if max(abs(d)) < 0.5
        s0 = pdf_MF_M2S_approx(d,0);
    elseif min(abs(d)) > 0.5
        s0 = pdf_MF_M2S_approx(d,1);
    else
        s0=[0;0;0];
    end
end

s=s0;
eps=1e-8;
alpha=0.05;
nf=2*eps;

NITER=0;
MAX_ITER=100;

sigma0=0.1;
sigma1=0.5;

while nf > eps && NITER < MAX_ITER
    
    NITER=NITER+1;

    f_stack=[];
    lambda_stack=[];

    if s(2)+s(3) > 1000
        [c,dc] = pdf_MF_normal(s,1,1);
        EQ = dc/c+1;
        f = EQ-d;
        df = get_J1(s);
    elseif s(1)+s(3) > 2000
        [c,dc] = pdf_MF_normal(s,1,1);
        EQ = dc/c+1;
        f = EQ-d;
        df = get_J2(s);
    else
        [c,dc,ddc] = pdf_MF_normal(s,1,1,1);
        EQ = dc/c+1;
        EQQ = zeros(3,3);
        for i = 1:3
            for j = 1:3
                EQQ(i,j) = (c+dc(i)+dc(j)+ddc(i,j))/c;
            end
        end

        df = EQQ-EQ*EQ';
        f = EQ-d;
    end
    
    s_dir=-df\f;
    f_stack=stack3(f_stack,f);
    lambda_stack=stack3(lambda_stack,0);
    
    lambda=1;
    s_trial=s+lambda*s_dir;
    while sum(s_trial<-1e-6) > 1
        lambda = lambda/2;
        s_trial=s+lambda*s_dir;
    end
    f_trial=func(s_trial,d,0,bool_scaled);
    f_stack=stack3(f_stack,f_trial);
    lambda_stack=stack3(lambda_stack,lambda);
    
    % polynomial line search: three-point parabolic method
    N_SUBITER=0;
    
    while norm(f_stack(end)) > (1-alpha*lambda)*norm(f) && N_SUBITER < MAX_ITER 
        
        N_SUBITER=N_SUBITER+1;
        if length(lambda_stack) < 3
            lambda=sigma1;
        else
            lam2=lambda_stack(2);
            lam3=lambda_stack(3);
            f2=f_stack(2);
            f3=f_stack(3);
            poly_coff_2=1/lam2/lam3/(lam2-lam3)*(lam3*(norm(f2)-norm(f))-lam2*(norm(f3)-norm(f)));
            if poly_coff_2 > 0
                poly_coff_1=1/(lam2-lam3)*(-lam3*(norm(f2)-norm(f))/lam2+lam2*(norm(f3)-norm(f))/lam3);
                lambda_t=-poly_coff_1/2/poly_coff_2;
                lambda=saturation(lambda_t,sigma0*lam3,sigma1*lam3)
            else
                lambda=sigma1*lam3;
            end

        end
        
        s_trial=s+lambda*s_dir;
        f_trial=func(s_trial,d,0,bool_scaled);
        f_stack=stack3(f_stack,f_trial);
        lambda_stack=stack3(lambda_stack,lambda);
   end
    
    s=s_trial;
    s=sort(s,'descend');
    nf=max(abs(f_trial));    
end

FVAL=f_trial;

if NITER == MAX_ITER || N_SUBITER == MAX_ITER
    disp('Warning: MAX iteration reached');
end

end

function Y=stack3(X,x)

Y=[X x];
n=size(Y,2);
if n > 3
    Y=Y(:,n-2:n);
end

end

function y=saturation(x,min_x,max_x)

if x <= min_x
    y=min_x;
elseif x >= max_x
    y=max_x;
else
    y=x;
end    
end


function varargout=func(s,d,bool_df,bool_scaled)
if ~bool_scaled
    if ~bool_df
        [c,dc]=pdf_MF_normal(s,0,1);
        f=dc-c*d;
        varargout{1}=f;
    else
        [c,dc,ddc]=pdf_MF_normal(s,0,1,1);
        f=dc-c*d;
        df=ddc-d*dc';
        varargout{1}=f;
        varargout{2}=df;
    end    
else
    if ~bool_df
        [c_bar,dc_bar]=pdf_MF_normal(s,1,1);
        f=dc_bar/c_bar-(d-ones(3,1));
        varargout{1}=f;
    else
        [c_bar,dc_bar,ddc_bar]=pdf_MF_normal(s,1,1,1);
        %f=dc_bar-c_bar*(d-ones(3,1));
        %df=ddc_bar-(d-ones(3,1))*dc_bar';        
        f=dc_bar/c_bar-(d-ones(3,1));
        df=ddc_bar/c_bar-dc_bar*dc_bar'/c_bar^2;
        varargout{1}=f;
        varargout{2}=df;
    end
end
end


function [ df ] = get_J1( s )

df(1,1) = 0.5*(1/(s(1)+s(2))^2+1/(s(1)+s(3))^2);
df(2,2) = 0.5*(1/(s(2)+s(1))^2+1/(s(2)+s(3))^2);
df(3,3) = 0.5*(1/(s(3)+s(1))^2+1/(s(3)+s(2))^2);
df(1,2) = 0.1/(s(1)+s(2))^2;
df(1,3) = 0.1/(s(1)+s(3))^2;
df(2,3) = 0.1/(s(2)+s(3))^2;
df(2,1) = df(1,2);
df(3,1) = df(1,3);
df(3,2) = df(2,3);

end


function [ df ] = get_J2( s )

I0 = besseli(0,s(2)+s(3),1);
I1 = besseli(1,s(2)+s(3),1);
I2 = besseli(2,s(2)+s(3),1);
I1I0 = I1/I0;
I2I0 = I2/I0;

s12 = 1/(s(1)+s(2));
s13 = 1/(s(1)+s(3));
s12s = s12^2;
s13s = s13^2;
s1213 = s12*s13;

df(1,1) = 0.5*(s12s+s13s);
df(2,2) = 0.5*(1+I2I0-2*I1I0^2)*(1-s12+0.5*s12s+1/12*s1213) +...
    1/4*(1+I2I0)*s12s + 1/8*(1-I2I0)*s1213;
df(3,3) = 0.5*(1+I2I0-2*I1I0^2)*(1-s13+0.5*s13s+1/12*s1213) +...
    1/4*(1+I2I0)*s13s + 1/8*(1-I2I0)*s1213;
df(1,2) = 0.5*I1I0*s12s;
df(1,3) = 0.5*I1I0*s13s;
df(2,3) = 0.5*(1+I2I0-2*I1I0^2)*(1-0.5*(s12+s13)+1/24*(3*s12s+3*s13s+2*s1213)) +...
    1/4*(I2I0-I1I0^2)*s1213;

df(2,1) = df(1,2);
df(3,1) = df(1,3);
df(3,2) = df(2,3);

end

