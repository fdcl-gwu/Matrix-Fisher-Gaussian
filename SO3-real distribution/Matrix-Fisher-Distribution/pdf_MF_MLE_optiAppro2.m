function [ F ] = pdf_MF_MLE_optiAppro2( ER, F0 )

% default initialization
if ~exist('F0','var')
    F0 = zeros(3);
end

F = F0;
[U,S,V] = psvd(F);
c = pdf_MF_normal(diag(S));
f = log(c)-trace(F*ER')

% step size
k = 1;

% iteration
i = 1;
while i == 1 || abs(f-fOld)>1e-8
    i = i+1;
    fOld = f;
    dc = pdf_MF_normal_deriv(diag(S));
    M = U*diag(dc/c)*V';
    dF = M-ER;
    F = F-k*dF;
    [U,S,V] = psvd(F);
    c = pdf_MF_normal(diag(S));
    f = log(c)-trace(F*ER')
end

end

