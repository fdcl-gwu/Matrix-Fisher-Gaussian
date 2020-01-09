function [ output ] = EQ( ind )

n = size(ind,1);

output = 0;
for k = 1:n
    comb = setC(n,k);
    Nc = length(comb);
    
    for nc = 1:Nc
        DcDs = {};
        DsDT = {};
        for nl = 1:k
            [DcDs,DsDT] = sumAlpha(DcDs,DsDT,ind(comb{nc}{nl},:));
        end
        
        output = output+sum(cellfun(@(x,y)x*y,DcDs,DsDT));
    end
end

end


function [ comb ] = setC( n, k )

if n==1 && k==1
    comb = {{1}};
    return;
elseif n==2 && k==1
    comb = {{[1,2]}};
    return;
end

if k>1
    comb = setC(n-1,k-1);
    Nc = length(comb);
    for nc = 1:Nc
        comb{nc} = [comb{nc},n];
    end
end

if n>k
    comb2 = setC(n-1,k);
    Nc = length(comb2);
    Nk = k;
    for nc = 1:Nc
        for nk = 1:Nk
            tempComb = comb2{nc};
            tempComb{nk} = [tempComb{nk},n];
            if ~exist('comb','var')
                comb = {tempComb};
            else
                comb = [comb,{tempComb}]; %#ok<AGROW>
            end
        end
    end
end

end


function [ DcDs, DsDT ] = sumAlpha( DcDs, DsDT, l )

if ~exist('DcDs','var') || isempty(DcDs)
    for nAlpha = 1:3
        DcDs{nAlpha} = sym(strcat('dcds',num2str(nAlpha)));
        DsDT{nAlpha} = dsdtT0(nAlpha,l);
    end
else
    N = length(DcDs);
    tempDcDs = cellfun(@char,DcDs,'UniformOutput',false);
    tempDsDT = DsDT;
    DcDs = cell(N*3,1);
    DsDT = cell(N*3,1);
    
    for n = 1:N
        for nAlpha = 1:3
            DcDs{3*(n-1)+nAlpha} = str2sym(strcat(tempDcDs{n},num2str(nAlpha)));
            DsDT{3*(n-1)+nAlpha} = tempDsDT{n}*dsdtT0(nAlpha,l);
        end
    end
end

end

