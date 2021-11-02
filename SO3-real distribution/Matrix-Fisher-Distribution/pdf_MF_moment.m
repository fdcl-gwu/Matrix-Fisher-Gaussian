function [ EQ, EQQ ] = pdf_MF_moment( s, bool_EQQ, bool_approx, s_threshold )

% if s are approximately equal, average them to eliminate numerical errors
% when calculating 1/(s(i)^2-s(j)^2)

e12 = abs(s(1)-s(2)) <= 1e-12*max([1,abs(s(1)),abs(s(2))]);
e23 = abs(s(2)-s(3)) <= 1e-12*max([1,abs(s(2)),abs(s(3))]);
if e12 && e23
    s(1:3) = mean(s);
elseif e12
    s(1:2) = mean(s(1:2));
elseif e23
    s(2:3) = mean(s(2:3));
end

% determine whether using highly concentrated approximations
if nargin < 3 || isempty(bool_approx)
    bool_approx = false;
end
if nargin < 4 || isempty(s_threshold)
    s_threshold = 1000;
end

if bool_approx && min([s(2)+s(3),s(1)+s(3),s(1)+s(2)])<s_threshold
    bool_approx = false;
end

% normalizing constant and derivatives
if nargin < 2 || isempty(bool_EQQ)
    bool_EQQ = false;
end

if ~bool_approx
    if ~bool_EQQ
        [c,dc] = pdf_MF_normal(s,1,1);
    else
        [c,dc,ddc] = pdf_MF_normal(s,1,1,1);
    end
else
    if ~bool_EQQ
        [c,dc] = pdf_MF_normal_approx(s,1,1);
    else
        [c,dc,ddc] = pdf_MF_normal_approx(s,1,1,1);
    end
end

%% first order moments
EQ = dc/c+1;
EQ = diag(EQ);

if ~bool_EQQ
    return;
end

%% second order moments
EQQ = zeros(9,9);

% EQiiQjj
for i = 1:3
    for j = 1:3
        ind1 = 3*(i-1)+i;
        ind2 = 3*(j-1)+j;
        EQQ(ind1,ind2) = (c+dc(i)+dc(j)+ddc(i,j))/c;
    end
end

% EQijQij
for i = 1:3
    for j = setdiff(1:3,i)
        ind1 = 3*(i-1)+j;
        if s(i)==s(j)
            EQQ(ind1,ind1) = 1/c/(2*s(i))*(c+dc(i)+s(i)*(ddc(i,i)-ddc(i,j)));
        else
            EQQ(ind1,ind1) = (1+dc(i)/c)*s(i)/(s(i)^2-s(j)^2)...
                -(1+dc(j)/c)*s(j)/(s(i)^2-s(j)^2);
        end
    end
end

% EQijQji
for i = 1:3
    for j = setdiff(1:3,i)
        ind1 = 3*(i-1)+j;
        ind2 = 3*(j-1)+i;
        if s(i)==s(j)
            EQQ(ind1,ind2) = 1/c/(2*s(i))*(-c-dc(i)+s(i)*(ddc(i,i)-ddc(i,j)));
        else
            EQQ(ind1,ind2) = (1+dc(i)/c)*s(j)/(s(i)^2-s(j)^2)...
                -(1+dc(j)/c)*s(i)/(s(i)^2-s(j)^2);
        end
    end
end

end

