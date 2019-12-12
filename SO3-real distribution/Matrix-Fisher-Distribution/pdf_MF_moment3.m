function [ EQQQ ] = pdf_MF_moment3( s, index, c, dc, ddc, dddc, scaled )
% index: if n-by-3, only calculate the third order moments specified by
% each row of index, return a n-by-1 vector; if empty, calculate all third
% order moments, return a 9-by-9-by-9 array
% c: if empty, calculate the normalising constant
% dc, ddc, dddc: if one of them is empty, calculate the derivatives
% scaled: default to true, if true, interpret c, dc, ddc, dddc as scaled
% versions.

if ~exist('scaled','var') || isempty(scaled)
    scaled = true;
end

if ~exist('c','var') || isempty(c)
    c = pdf_MF_normal(s,scaled);
end
if ~exist('dc','var') || isempty(dc) ||...
        ~exist('ddc','var') || isempty(ddc) ||...
        ~exist('dddc','var') || isempty(dddc)
    [dc,ddc,dddc] = pdf_MF_normal_deriv(s,true,scaled);
end

% calculate dc/c, ddc/c, dddc/c
if scaled
    % DO NOT CHANGE THE ORDER!
    % dddc/c
    for i = 1:3
        for j = 1:3
            for k = 1:3
                dddc(i,j,k) = (dddc(i,j,k)+ddc(i,j)+ddc(i,k)+ddc(j,k)...
                    +dc(i)+dc(j)+dc(k)+c)/c;
            end
        end
    end
    
    % ddc/c
    for i = 1:3
        for j = 1:3
            ddc(i,j) = (ddc(i,j)+dc(i)+dc(j)+c)/c;
        end
    end
    
    % dc/c
    for i = 1:3
        dc(i) = (dc(i)+c)/c;
    end
else
    dc = dc/c;
    ddc = ddc/c;
    dddc = dddc/c;
end

if ~exist('index','var') || isempty(index)
    EQQQ = zeros(9,9,9);
    for i = 1:3
        for j = 1:3
            for k = 1:3
                % iijjkk
                EQQQ(3*(i-1)+i,3*(j-1)+j,3*(k-1)+k) = dddc(i,j,k);

                if j~=k
                    if i~=j && i~=k
                        % iijkjk
                        EQQQ(3*(i-1)+i,3*(j-1)+k,3*(j-1)+k) = ...
                            (ddc(i,k)*s(k)-ddc(i,j)*s(j))/(s(k)^2-s(j)^2);
                        % iijkkj
                        EQQQ(3*(i-1)+i,3*(j-1)+k,3*(k-1)+j) = ...
                            (ddc(i,k)*s(j)-ddc(i,j)*s(k))/(s(k)^2-s(j)^2);

                        % jkiijk & jkiikj
                        EQQQ(3*(j-1)+k,3*(i-1)+i,3*(j-1)+k) = EQQQ(3*(i-1)+i,3*(j-1)+k,3*(j-1)+k);
                        EQQQ(3*(j-1)+k,3*(i-1)+i,3*(k-1)+j) = EQQQ(3*(i-1)+i,3*(j-1)+k,3*(k-1)+j);

                        % jkjkii & jkkjii
                        EQQQ(3*(j-1)+k,3*(j-1)+k,3*(i-1)+i) = EQQQ(3*(i-1)+i,3*(j-1)+k,3*(j-1)+k);
                        EQQQ(3*(j-1)+k,3*(k-1)+j,3*(i-1)+i) = EQQQ(3*(i-1)+i,3*(j-1)+k,3*(k-1)+j);
                    elseif i==j
                        % jjjkjk
                        EQQQ(3*(j-1)+j,3*(j-1)+k,3*(j-1)+k) = ...
                            ((ddc(j,k)*s(k)-ddc(j,j)*s(j))/(s(k)^2-s(j)^2)+...
                            (-dc(j)*(s(j)^2+s(k)^2)+dc(k)*2*s(j)*s(k))/(s(j)^2-s(k)^2)^2);
                        % jjjkkj
                        EQQQ(3*(j-1)+j,3*(j-1)+k,3*(k-1)+j) = ...
                            ((ddc(j,k)*s(j)-ddc(j,j)*s(k))/(s(k)^2-s(j)^2)+...
                            (-dc(j)*2*s(j)*s(k)+dc(k)*(s(j)^2+s(k)^2))/(s(j)^2-s(k)^2)^2);

                        % jkjjjk & jkjjkj
                        EQQQ(3*(j-1)+k,3*(j-1)+j,3*(j-1)+k) = EQQQ(3*(j-1)+j,3*(j-1)+k,3*(j-1)+k);
                        EQQQ(3*(j-1)+k,3*(j-1)+j,3*(k-1)+j) = EQQQ(3*(j-1)+j,3*(j-1)+k,3*(k-1)+j);

                        % jkjkjj & jkkjjj
                        EQQQ(3*(j-1)+k,3*(j-1)+k,3*(j-1)+j) = EQQQ(3*(j-1)+j,3*(j-1)+k,3*(j-1)+k);
                        EQQQ(3*(j-1)+k,3*(k-1)+j,3*(j-1)+j) = EQQQ(3*(j-1)+j,3*(j-1)+k,3*(k-1)+j);
                    else
                        % kkjkjk & kkjkkj
                        EQQQ(3*(k-1)+k,3*(j-1)+k,3*(j-1)+k) = ...
                            ((ddc(k,j)*s(j)-ddc(k,k)*s(k))/(s(j)^2-s(k)^2)+...
                            (-dc(k)*(s(k)^2+s(j)^2)+dc(j)*2*s(k)*s(j))/(s(k)^2-s(j)^2)^2);
                        EQQQ(3*(k-1)+k,3*(j-1)+k,3*(k-1)+j) = ...
                            ((ddc(k,j)*s(k)-ddc(k,k)*s(j))/(s(j)^2-s(k)^2)+...
                            (-dc(k)*2*s(k)*s(j)+dc(j)*(s(k)^2+s(j)^2))/(s(k)^2-s(j)^2)^2);

                        % jkkkjk & jkkkkj
                        EQQQ(3*(j-1)+k,3*(k-1)+k,3*(j-1)+k) = EQQQ(3*(k-1)+k,3*(j-1)+k,3*(j-1)+k);
                        EQQQ(3*(j-1)+k,3*(k-1)+k,3*(k-1)+j) = EQQQ(3*(k-1)+k,3*(j-1)+k,3*(k-1)+j);

                        % jkjkkk & jkkjkk
                        EQQQ(3*(j-1)+k,3*(j-1)+k,3*(k-1)+k) = EQQQ(3*(k-1)+k,3*(j-1)+k,3*(j-1)+k);
                        EQQQ(3*(j-1)+k,3*(k-1)+j,3*(k-1)+k) = EQQQ(3*(k-1)+k,3*(j-1)+k,3*(k-1)+j);
                    end
                end

                if i~=j && i~=k && j~=k
                    % ijjkki
                    EQQQ(3*(i-1)+j,3*(j-1)+k,3*(k-1)+i) = dc(i)*s(j)*s(k)/(s(i)^2-s(j)^2)/(s(i)^2-s(k)^2)...
                        + dc(j)*s(i)*s(k)/(s(j)^2-s(i)^2)/(s(j)^2-s(k)^2) + dc(k)*s(i)*s(j)/(s(k)^2-s(i)^2)/(s(k)^2-s(j)^2);
                    % ijjkik
                    EQQQ(3*(i-1)+j,3*(j-1)+k,3*(i-1)+k) = dc(i)*s(i)*s(j)/(s(j)^2-s(i)^2)/(s(k)^2-s(i)^2)...
                        + dc(j)*s(j)*s(j)/(s(i)^2-s(j)^2)/(s(k)^2-s(j)^2) + dc(k)*s(j)*s(k)/(s(j)^2-s(k)^2)/(s(i)^2-s(k)^2);

                    % ijkijk & ijikjk
                    EQQQ(3*(i-1)+j,3*(k-1)+i,3*(j-1)+k) = EQQQ(3*(i-1)+j,3*(j-1)+k,3*(k-1)+i);
                    EQQQ(3*(i-1)+j,3*(i-1)+k,3*(j-1)+k) = EQQQ(3*(i-1)+j,3*(j-1)+k,3*(i-1)+k);

                    % jkijki & jkijik
                    EQQQ(3*(j-1)+k,3*(i-1)+j,3*(k-1)+i) = EQQQ(3*(i-1)+j,3*(j-1)+k,3*(k-1)+i);
                    EQQQ(3*(j-1)+k,3*(i-1)+j,3*(i-1)+k) = EQQQ(3*(i-1)+j,3*(j-1)+k,3*(i-1)+k);
                    
                    % jkkiij & jkikij
                    EQQQ(3*(j-1)+k,3*(k-1)+i,3*(i-1)+j) = EQQQ(3*(i-1)+j,3*(j-1)+k,3*(k-1)+i);
                    EQQQ(3*(j-1)+k,3*(i-1)+k,3*(i-1)+j) = EQQQ(3*(i-1)+j,3*(j-1)+k,3*(i-1)+k);

                    % kiijjk & ikijjk
                    EQQQ(3*(k-1)+i,3*(i-1)+j,3*(j-1)+k) = EQQQ(3*(i-1)+j,3*(j-1)+k,3*(k-1)+i);
                    EQQQ(3*(i-1)+k,3*(i-1)+j,3*(j-1)+k) = EQQQ(3*(i-1)+j,3*(j-1)+k,3*(i-1)+k);
                    
                    % kijkij & ikjkij
                    EQQQ(3*(k-1)+i,3*(j-1)+k,3*(i-1)+j) = EQQQ(3*(i-1)+j,3*(j-1)+k,3*(k-1)+i);
                    EQQQ(3*(i-1)+k,3*(j-1)+k,3*(i-1)+j) = EQQQ(3*(i-1)+j,3*(j-1)+k,3*(i-1)+k);
                end
             end
        end
    end
else
    N = size(index,1);
    EQQQ = zeros(N,1);
    
    index1 = index(:,[1,3,5]);
    index2 = index(:,[2,4,6]);
    for n = 1:N
        if rem(sum(1==index(n,:)),2)~=0 ||...
                rem(sum(2==index(n,:)),2)~=0 ||...
                rem(sum(3==index(n,:)),2)~=0
            % odd repeatation index
            EQQQ(n) = 0;
        elseif ~sum(index1~=index2)
            % iijjkk
            EQQQ(n) = dddc(index1(n,1),index1(n,2),index1(n,3));
        else
            repeat1 = find(index1==index2);
            if isempty(repeat1)
                % there isn't ii
                i = index1(1);
                j = index2(1);
                k = setdiff([1,2,3],[i,j]);
                if length(unique(index1))==3
                    % ijjkki
                    EQQQ(n) = dc(i)*s(j)*s(k)/(s(i)^2-s(j)^2)/(s(i)^2-s(k)^2)...
                        + dc(j)*s(i)*s(k)/(s(j)^2-s(i)^2)/(s(j)^2-s(k)^2) + dc(k)*s(i)*s(j)/(s(k)^2-s(i)^2)/(s(k)^2-s(j)^2);
                elseif intersect(unique(index1),unique(index2))==i
                    % jiikjk
                    EQQQ(n) = dc(i)*s(i)*s(i)/(s(j)^2-s(i)^2)/(s(k)^2-s(i)^2)...
                        + dc(j)*s(i)*s(j)/(s(i)^2-s(j)^2)/(s(k)^2-s(j)^2) + dc(k)*s(i)*s(k)/(s(i)^2-s(k)^2)/(s(j)^2-s(k)^2);
                elseif intersect(unique(index1),unique(index2))==j
                    % ijjkik
                    EQQQ(n) = dc(i)*s(i)*s(j)/(s(j)^2-s(i)^2)/(s(k)^2-s(i)^2)...
                        + dc(j)*s(j)*s(j)/(s(i)^2-s(j)^2)/(s(k)^2-s(j)^2) + dc(k)*s(j)*s(k)/(s(j)^2-s(k)^2)/(s(i)^2-s(k)^2);
                else
                    % ikkjij
                    EQQQ(n) = dc(i)*s(i)*s(k)/(s(j)^2-s(i)^2)/(s(k)^2-s(i)^2)...
                        + dc(j)*s(j)*s(k)/(s(i)^2-s(j)^2)/(s(k)^2-s(j)^2) + dc(k)*s(k)*s(k)/(s(j)^2-s(k)^2)/(s(i)^2-s(k)^2);
                end
            else
                % there is ii
                norepeat = setdiff([1,2,3],repeat1);
                i = index1(repeat1);
                j = index1(norepeat(1));
                k = index2(norepeat(1));
                if index1(norepeat(1))==index1(norepeat(2))
                    if i==j
                        % jjjkjk
                        EQQQ(n) = (ddc(j,k)*s(k)-ddc(j,j)*s(j))/(s(k)^2-s(j)^2)+...
                            (-dc(j)*(s(j)^2+s(k)^2)+dc(k)*2*s(j)*s(k))/(s(k)^2-s(j)^2)^2;
                    elseif i==k
                        % kkjkjk
                        EQQQ(n) = (ddc(k,j)*s(j)-ddc(k,k)*s(k))/(s(j)^2-s(k)^2)+...
                            (-dc(k)*(s(k)^2+s(j)^2)+dc(j)*2*s(k)*s(j))/(s(j)^2-s(k)^2)^2;
                    else
                        % iijkjk
                        EQQQ(n) = (ddc(i,k)*s(k)-ddc(i,j)*s(j))/(s(k)^2-s(j)^2);
                    end
                else
                    if i==j
                        % jjjkkj
                        EQQQ(n) = (ddc(j,k)*s(j)-ddc(j,j)*s(k))/(s(k)^2-s(j)^2)+...
                            (-dc(j)*2*s(j)*s(k)+dc(k)*(s(j)^2+s(k)^2))/(s(k)^2-s(j)^2)^2;
                    elseif i==k
                        % kkjkkj
                        EQQQ(n) = (ddc(k,j)*s(k)-ddc(k,k)*s(j))/(s(j)^2-s(k)^2)+...
                            (-dc(k)*2*s(k)*s(j)+dc(j)*(s(k)^2+s(j)^2))/(s(j)^2-s(k)^2)^2;
                    else
                        % iijkkj
                        EQQQ(n) = (ddc(i,k)*s(j)-ddc(i,j)*s(k))/(s(k)^2-s(j)^2);
                    end
                end
            end
        end
    end
end

end

