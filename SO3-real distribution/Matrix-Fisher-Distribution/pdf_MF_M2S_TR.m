function [ s, n_iter ] = pdf_MF_M2S_TR( d, s0, options )

if nargin < 2 || isempty(s0)
    if max(abs(d)) < 0.5
        s0 = pdf_MF_M2S_approx(d,0);
    elseif min(abs(d)) > 0.5
        s0 = pdf_MF_M2S_approx(d,1);
    else
        s0=[0;0;0];
    end
end

% lsqlin options
if nargin < 3 || isempty(options)
    options = optimoptions('lsqlin','display','off','algorithm','active-set');
end

% trust region parameters
eta_v = 0.9;
eta_s = 0;
gamma_i = 2;
gamma_d = 0.5;

Delta = [100;100;100];
eps = 1e-16;

iter_max = 100;

% initial step
s = s0;

[c,dc,ddc]=pdf_MF_normal(s,1,1,1);
EQ = dc/c+1;
f = sum((EQ-d).^2);

n_iter = 0;

while f > eps
    % calculate jacobias
    if s(2)+s(3) > 1000   
        J = get_J1(s);
    elseif s(1)+s(3) > 1000
        J = get_J2(s);
    else
        EQQ = zeros(3,3);
        for i = 1:3
            for j = 1:3
                EQQ(i,j) = (c+dc(i)+dc(j)+ddc(i,j))/c;
            end
        end

        J = EQQ-EQ*EQ';
    end
    
    % solve constrained LS
    dg = EQ-d;
    ds = J\(-dg);
    if sum(abs(ds)<Delta) < 3
        [ds,f_LS] = lsqlin(J,-dg,[],[],[],[],-Delta,Delta,ds,options);
    else
        f_LS = sum((J*ds+dg).^2);
    end
    
    % evaluate trust region
    s_new = sort(s+ds,'descend');
    [c_new,dc_new,ddc_new]=pdf_MF_normal(s_new,1,1,1);
    EQ_new = dc_new/c_new+1;
    f_new = sum((EQ_new-d).^2);
    
    rho = (f-f_new)/(f-f_LS);
    
    % update
    if rho >= eta_v
        Delta = gamma_i*Delta;
    elseif rho < eta_s
        Delta = gamma_d*Delta;
    end
    
    if rho >= eta_s
        s = s_new;
        c = c_new;
        dc = dc_new;
        ddc = ddc_new;
        EQ = EQ_new;
        f = f_new;
    end
    
    n_iter = n_iter+1;
    if n_iter > iter_max
        warning('warning: M2S maximum iteration reached');
        break;
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
    1/4*(1-I1I0^2)*s1213;

df(2,1) = df(1,2);
df(3,1) = df(1,3);
df(3,2) = df(2,3);

end

