function [ EQQ ] = pdf_MF_moment2( s )

c = pdf_MF_normal(s,true);
[dc,ddc] = pdf_MF_normal_deriv(s,true,true);

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

