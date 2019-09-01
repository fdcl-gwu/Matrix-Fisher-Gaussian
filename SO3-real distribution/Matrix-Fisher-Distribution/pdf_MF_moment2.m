function [ EQQ ] = pdf_MF_moment2( S )

c = pdf_MF_normal(S,true);
[dc,ddc] = pdf_MF_normal_deriv(S,true,true);

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
        EQQ(ind1,ind1) = (1+dc(i)/c)*S(i)/(S(i)^2-S(j)^2)...
            -(1+dc(j)/c)*S(j)/(S(i)^2-S(j)^2);
    end
end

% EQijQji
for i = 1:3
    for j = setdiff(1:3,i)
        ind1 = 3*(i-1)+j;
        ind2 = 3*(j-1)+i;
        EQQ(ind1,ind2) = (1+dc(i)/c)*S(j)/(S(i)^2-S(j)^2)...
            -(1+dc(j)/c)*S(i)/(S(i)^2-S(j)^2);
    end
end

end

