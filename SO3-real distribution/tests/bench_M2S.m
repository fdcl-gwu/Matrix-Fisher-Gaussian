clear;
close all;

addpath('../Matrix-Fisher-Distribution');

% test cases
Ns = 20;

s = zeros(3,Ns,Ns,2*Ns);
s1 = logspace(-1,5,Ns);

for i = 1:Ns
    s(1,i,:,:) = s1(i);
    s2 = logspace(-log10(s1(i))-6,log10(s1(i)),Ns);
    
    for j = 1:Ns
        s(2,i,j,:) = s2(j);
        s3 = [-flip(logspace(log10(s2(j))-6,log10(s2(j)),Ns)),...
            logspace(log10(s2(j))-6,log10(s2(j)),Ns)];
        s(3,i,j,:) = permute(s3,[1,3,4,2]);
    end
end

% calculate
sinv_LS = zeros(3,Ns,Ns,2*Ns);
error_LS = zeros(Ns,Ns,2*Ns);
NIter_LS = zeros(Ns,Ns,2*Ns);
success_LS = true(Ns,Ns,2*Ns);
time_LS = zeros(Ns,Ns,2*Ns);

sinv_TR = zeros(3,Ns,Ns,2*Ns);
error_TR = zeros(Ns,Ns,2*Ns);
NIter_TR = zeros(Ns,Ns,2*Ns);
success_TR = true(Ns,Ns,2*Ns);
time_TR = zeros(Ns,Ns,2*Ns);

options = optimoptions('lsqlin','display','off');

for i = 15:Ns
    for j = 1:Ns
        for k = 1:2*Ns
            d = pdf_MF_moment(s(:,i,j,k));
            
            if i==3 && j==19 && k==1
                a = 1;
            end
            
            try
                tic;
                [sinv_LS(:,i,j,k),~,NIter_LS(i,j,k)] = pdf_MF_M2S(d);
                time_LS(i,j,k) = toc;
                error_LS(i,j,k) = max(abs(s(:,i,j,k)-sinv_LS(:,i,j,k)))/...
                    sqrt(sum(s(:,i,j,k).^2));
                if sum(isfinite(sinv_LS(:,i,j,k))) < 3
                    success_LS(i,j,k) = false;
                end
            catch
                success_LS(i,j,k) = false;
            end
            
%             try
%                 tic;
%                 [sinv_TR(:,i,j,k),NIter_TR(i,j,k)] = pdf_MF_M2S_TR(d,[],options);
%                 time_TR(i,j,k) = toc;
%                 error_TR(i,j,k) = max(abs(s(:,i,j,k)-sinv_TR(:,i,j,k)))/...
%                     sqrt(sum(s(:,i,j,k).^2));
%                 if sum(isfinite(sinv_TR(:,i,j,k))) < 3
%                     success_TR(i,j,k) = false;
%                 end
%             catch
%                 success_TR(i,j,k) = false;
%             end
        end
    end
end

% save data
path = 'D:\result-SO3Euclid\bench_M2S\10-28-2021';
if ~exist(path,'dir')
    mkdir(path);
end

filename = strcat(path,'\data.mat');
save(filename,'s','sinv_LS','error_LS','NIter_LS','success_LS','time_LS',...
    'sinv_TR','error_TR','NIter_TR','success_TR','time_TR');

rmpath('../Matrix-Fisher-Distribution');


