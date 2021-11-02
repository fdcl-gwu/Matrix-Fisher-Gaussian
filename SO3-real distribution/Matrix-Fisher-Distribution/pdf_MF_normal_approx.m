function [c_return,dc_return,ddc_return] = pdf_MF_normal_approx(s,bool_scaled,bool_dc,bool_ddc)

% if bool_scaled is not defined, then set it false
if nargin < 2
    bool_scaled=false;
end
if nargin < 3
    bool_dc=false;
end
if nargin < 4
    bool_ddc=false;
end
if bool_ddc
    bool_dc = true;
end

%% normalizing constant
if ~bool_scaled
    % return the normalizing constant without any scaling
    c = exp(sum(s))/(8*pi*(s(1)+s(2))*(s(2)+s(3))*(s(1)+s(3)));
    c_return = c;
else
    % return the normalizing constant scaled by exp(-sum(s))
    c_bar = 1/(8*pi*(s(1)+s(2))*(s(2)+s(3))*(s(1)+s(3)));
    c_return = c_bar;
end

if ~bool_dc
    return;
end

%% first order derivative
if ~bool_scaled
    % derivatives of the normalizing constant
    dc(1) = c*(1-0.5*(1/(s(1)+s(2))+1/(s(1)+s(3))));
    dc(2) = c*(1-0.5*(1/(s(2)+s(1))+1/(s(2)+s(3))));
    dc(3) = c*(1-0.5*(1/(s(3)+s(1))+1/(s(3)+s(2))));
    dc_return = dc;
else
    % derivatives of the scaled normalizing constant
    dc_bar(1) = -0.5*c_bar*(1/(s(1)+s(2))+1/(s(1)+s(3)));
    dc_bar(2) = -0.5*c_bar*(1/(s(2)+s(1))+1/(s(2)+s(3)));
    dc_bar(3) = -0.5*c_bar*(1/(s(3)+s(1))+1/(s(3)+s(2)));
    dc_return = dc_bar;
end

if ~bool_ddc
    return;
end

%% second order derivative
if ~bool_scaled
    A = zeros(9,9);
    b = zeros(9,1);

    for i = 1:3
        for j = 1:3
            k = setdiff(1:3,[i,j]);
            if i==j
                if abs(s(i))~=abs(s(k(1))) && abs(s(i))~=abs(s(k(2)))
                    A(3*(i-1)+j,3*(i-1)+j) = 1;
                    b(3*(i-1)+j) = c-(-dc(i)*s(i)+dc(k(1))*s(k(1)))/(s(k(1))^2-s(i)^2)-(-dc(i)*s(i)+dc(k(2))*s(k(2)))/(s(k(2))^2-s(i)^2);
                elseif abs(s(i))~=abs(s(k(1))) && abs(s(i))==abs(s(k(2))) && s(i)~=0
                    A(3*(i-1)+j,3*(i-1)+j) = 3/2;
                    A(3*(i-1)+j,3*(i-1)+k(2)) = -1/2*sign(s(i)*s(k(2)));
                    b(3*(i-1)+j) = c-(-dc(i)*s(i)+dc(k(1))*s(k(1)))/(s(k(1))^2-s(i)^2)-dc(i)/2/s(i);
                elseif abs(s(i))~=abs(s(k(1))) && s(i)==s(k(2)) && s(i)==0
                    A(3*(i-1)+j,3*(i-1)+j) = 2;
                    b(3*(i-1)+j) = c-(-dc(i)*s(i)+dc(k(1))*s(k(1)))/(s(k(1))^2-s(i)^2);
                elseif abs(s(i))==abs(s(k(1))) && abs(s(i))~=abs(s(k(2))) && s(i)~=0
                    A(3*(i-1)+j,3*(i-1)+j) = 3/2;
                    A(3*(i-1)+j,3*(i-1)+k(1)) = -1/2*sign(s(i)*s(k(1)));
                    b(3*(i-1)+j) = c-dc(i)/2/s(i)-(-dc(i)*s(i)+dc(k(2))*s(k(2)))/(s(k(2))^2-s(i)^2);
                elseif abs(s(i))==abs(s(k(1))) && abs(s(i))==abs(s(k(2))) && s(i)~=0
                    A(3*(i-1)+j,3*(i-1)+j) = 2;
                    A(3*(i-1)+j,3*(i-1)+k(1)) = -1/2*sign(s(i)*s(k(1)));
                    A(3*(i-1)+j,3*(i-1)+k(2)) = -1/2*sign(s(i)*s(k(2)));
                    b(3*(i-1)+j) = c-dc(i)/s(i);
                elseif s(i)==s(k(1)) && abs(s(i))~=abs(s(k(2))) && s(i)==0
                    A(3*(i-1)+j,3*(i-1)+j) = 2;
                    b(3*(i-1)+j) = c-(-dc(i)*s(i)+dc(k(2))*s(k(2)))/(s(k(2))^2-s(i)^2);
                else
                    A(3*(i-1)+j,3*(i-1)+j) = 3;
                    b(3*(i-1)+j)=c;
                end
            else
                if abs(s(i))~=abs(s(j))
                    A(3*(i-1)+j,3*(i-1)+j) = 1;
                    b(3*(i-1)+j) = dc(k)+(-dc(i)*s(j)+dc(j)*s(i))/(s(j)^2-s(i)^2);
                elseif abs(s(i))==abs(s(j)) && s(i)~=0
                    A(3*(i-1)+j,3*(i-1)+j) = 3/2;
                    A(3*(i-1)+j,3*(i-1)+i) = -1/2*sign(s(i)*s(j));
                    b(3*(i-1)+j) = dc(k)-dc(j)/2/s(i);
                else
                    A(3*(i-1)+j,3*(i-1)+j) = 2;
                    b(3*(i-1)+j) = dc(k);
                end
            end
        end
    end

    ddc = A\b;
    ddc = [ddc(1:3),ddc(4:6),ddc(7:9)];
    ddc_return = ddc;
else
    A = zeros(9,9);
    b = zeros(9,1);

    for i = 1:3
        for j = 1:3
            k = setdiff(1:3,[i,j]);
            if i==j
                if abs(s(i))~=abs(s(k(1))) && abs(s(i))~=abs(s(k(2)))
                    A(3*(i-1)+j,3*(i-1)+j) = 1;
                    b(3*(i-1)+j) = -2*dc_bar(i) - c_bar/(s(i)+s(k(1)))+dc_bar(i)*s(i)/(s(k(1))^2-s(i)^2)-dc_bar(k(1))*s(k(1))/(s(k(1))^2-s(i)^2)...
                        - c_bar/(s(i)+s(k(2)))+dc_bar(i)*s(i)/(s(k(2))^2-s(i)^2)-dc_bar(k(2))*s(k(2))/(s(k(2))^2-s(i)^2);
                elseif abs(s(i))~=abs(s(k(1))) && abs(s(i))==abs(s(k(2))) && s(i)~=0
                    sig = sign(s(i)*s(k(2)));
                    A(3*(i-1)+j,3*(i-1)+j) = 3/2;
                    A(3*(i-1)+j,3*(i-1)+k(2)) = -1/2*sig;
                    b(3*(i-1)+j) = -2*dc_bar(i) - c_bar/(s(i)+s(k(1)))+dc_bar(i)*s(i)/(s(k(1))^2-s(i)^2)-dc_bar(k(1))*s(k(1))/(s(k(1))^2-s(i)^2)...
                        - (1/2-sig/2+1/2/s(i))*c_bar-(1/2/s(i)+1-sig/2)*dc_bar(i)+sig*dc_bar(k(2))/2;
                elseif abs(s(i))~=abs(s(k(1))) && s(i)==s(k(2)) && s(i)==0
                    A(3*(i-1)+j,3*(i-1)+j) = 2;
                    b(3*(i-1)+j) = -2*dc_bar(i) - c_bar/(s(i)+s(k(1)))+dc_bar(i)*s(i)/(s(k(1))^2-s(i)^2)-dc_bar(k(1))*s(k(1))/(s(k(1))^2-s(i)^2)...
                        - c_bar-2*dc_bar(i);
                elseif abs(s(i))==abs(s(k(1))) && abs(s(i))~=abs(s(k(2))) && s(i)~=0
                    sig = sign(s(i)*s(k(1)));
                    A(3*(i-1)+j,3*(i-1)+j) = 3/2;
                    A(3*(i-1)+j,3*(i-1)+k(1)) = -1/2*sign(s(i)*s(k(1)));
                    b(3*(i-1)+j) = -2*dc_bar(i) - (1/2-sig/2+1/2/s(i))*c_bar-(1/2/s(i)+1-sig/2)*dc_bar(i)+sig*dc_bar(k(1))/2 ...
                        - c_bar/(s(i)+s(k(2)))+dc_bar(i)*s(i)/(s(k(2))^2-s(i)^2)-dc_bar(k(2))*s(k(2))/(s(k(2))^2-s(i)^2);
                elseif abs(s(i))==abs(s(k(1))) && abs(s(i))==abs(s(k(2))) && s(i)~=0
                    sig1 = sign(s(i)*s(k(1)));
                    sig2 = sign(s(i)*s(k(2)));
                    A(3*(i-1)+j,3*(i-1)+j) = 2;
                    A(3*(i-1)+j,3*(i-1)+k(1)) = -1/2*sig1;
                    A(3*(i-1)+j,3*(i-1)+k(2)) = -1/2*sig2;
                    b(3*(i-1)+j) = -2*dc_bar(i) - (1/2-sig1/2+1/2/s(i))*c_bar-(1/2/s(i)+1-sig1/2)*dc_bar(i)+sig1*dc_bar(k(1))/2 ...
                        - (1/2-sig2/2+1/2/s(i))*c_bar-(1/2/s(i)+1-sig2/2)*dc_bar(i)+sig2*dc_bar(k(2))/2;
                elseif s(i)==s(k(1)) && abs(s(i))~=abs(s(k(2))) && s(i)==0
                    A(3*(i-1)+j,3*(i-1)+j) = 2;
                    b(3*(i-1)+j) = -2*dc_bar(i) - c_bar-2*dc_bar(i)...
                        - c_bar/(s(i)+s(k(2)))+dc_bar(i)*s(i)/(s(k(2))^2-s(i)^2)-dc_bar(k(2))*s(k(2))/(s(k(2))^2-s(i)^2);
                else
                    A(3*(i-1)+j,3*(i-1)+j) = 3;
                    b(3*(i-1)+j) = -2*dc_bar(i) - c_bar-2*dc_bar(i) - c_bar-2*dc_bar(i);
                end
            else
                if abs(s(i))~=abs(s(j))
                    A(3*(i-1)+j,3*(i-1)+j) = 1;
                    b(3*(i-1)+j) = -c_bar/(s(i)+s(j)) - (1+s(j)/(s(j)^2-s(i)^2))*dc_bar(i)...
                        - (1-s(i)/(s(j)^2-s(i)^2))*dc_bar(j) + dc_bar(k);
                elseif abs(s(i))==abs(s(j)) && s(i)~=0
                    sig = sign(s(i)*s(j));
                    A(3*(i-1)+j,3*(i-1)+j) = 3/2;
                    A(3*(i-1)+j,3*(i-1)+i) = -1/2*sig;
                    b(3*(i-1)+j) = -(1/2-sig/2+1/2/s(i))*c_bar - (3/2-sig)*dc_bar(i) - (3/2+1/2/s(i))*dc_bar(j) + dc_bar(k);
                else
                    A(3*(i-1)+j,3*(i-1)+j) = 2;
                    b(3*(i-1)+j) = -c_bar - 2*dc_bar(i) - 2*dc_bar(j) + dc_bar(k);
                end
            end
        end
    end

    ddc_bar = A\b;
    ddc_bar = [ddc_bar(1:3),ddc_bar(4:6),ddc_bar(7:9)];
    ddc_return = ddc_bar;
end

end

