function [ output, MatlabStr, texStr ] = simplifyQ( input )

input = simplify(input);
input = expand(input);

%% unique combinations
unique2 = {'1111','1122','1133','2222','2233','3333',...
    '1212','1221','1313','1331','2323','2332'};
unique2Sym = cell(12,1);
unique2Matlab = cell(12,1);
unique2Tex = cell(12,1);
for i = 1:12
    unique2Sym{i} = strcat('Q',unique2{i}([1,2]),'*Q',unique2{i}([3,4]));
    unique2Sym{i} = str2sym(unique2Sym{i});
    
    unique2Matlab{i} = strcat('EQQ(',num2str(3*(str2double(unique2{i}(1))-1)+str2double(unique2{i}(2))),...
        ',',num2str(3*(str2double(unique2{i}(3))-1)+str2double(unique2{i}(4))),')');
    
    if strcmp(unique2{i}(1:2),unique2{i}(3:4))
        unique2Tex{i} = strcat('\expect{Q_{',num2str(unique2{i}(1:2)),'}^2}');
    else
        unique2Tex{i} = strcat('\expect{Q_{',num2str(unique2{i}(1:2)),'}Q_{',...
            num2str(unique2{i}(3:4)),'}}');
    end
end

unique3 = {'111111','111122','111133','112222','112233','113333',...
    '222222','222233','223333','333333',...
    '111212','111221','111313','111331','112323','112332',...
    '221212','221221','221313','221331','222323','222332',...
    '331212','331221','331313','331331','332323','332332',...
    '122331','211323','122313','133212'};
unique3Sym = cell(32,1);
unique3Matlab = cell(12,1);
unique3Tex = cell(12,1);
for i = 1:32
    unique3Sym{i} = strcat('Q',unique3{i}([1,2]),'*Q',unique3{i}([3,4]),...
        '*Q',unique3{i}([5,6]));
    unique3Sym{i} = str2sym(unique3Sym{i});
    
    unique3Matlab{i} = strcat('EQQQ(',num2str(3*(str2double(unique3{i}(1))-1)+str2double(unique3{i}(2))),...
        ',',num2str(3*(str2double(unique3{i}(3))-1)+str2double(unique3{i}(4))),...
        ',',num2str(3*(str2double(unique3{i}(5))-1)+str2double(unique3{i}(6))),')');
    
    if strcmp(unique3{i}(1:2),unique3{i}(3:4)) && strcmp(unique3{i}(1:2),unique3{i}(5:6))
        unique3Tex{i} = strcat('\expect{Q_{',num2str(unique3{i}(1:2)),'}^3}');
    elseif strcmp(unique3{i}(1:2),unique3{i}(3:4))
        unique3Tex{i} = strcat('\expect{Q_{',num2str(unique3{i}(1:2)),'}^2Q_{',...
            num2str(unique3{i}(5:6)),'}}');
    elseif strcmp(unique3{i}(3:4),unique3{i}(5:6))
        unique3Tex{i} = strcat('\expect{Q_{',num2str(unique3{i}(1:2)),'}Q_{',...
            num2str(unique3{i}(3:4)),'}^2}');
    else
        unique3Tex{i} = strcat('\expect{Q_{',num2str(unique3{i}(1:2)),'}Q_{',...
            num2str(unique3{i}(3:4)),'}Q_{',num2str(unique3{i}(5:6)),'}}');
    end
end

%%
[N1,N2] = size(input);
output = sym(zeros(N1,N2));
MatlabStr = cell(N1,N2);
MatlabStr(:) = {''};
texStr = cell(N1,N2);
texStr(:) = {''};
for n1 = 1:N1
    for n2 = 1:N2
        terms = children(input(n1,n2));
        
        %% delete zeros
        N = length(terms);
        isZero = false(N,1);
        subQ = cell(N,1);
        termStr = cell(N,1);
        for n = 1:N
            termStr{n} = char(terms{n});
            indQ = strfind(termStr{n},'Q');
            
            indN = length(indQ);
            subQ{n} = [];
            for indn = 1:indN
                if strcmp(termStr{n}(indQ(indn)+3),'^')
                    power = str2double(termStr{n}(indQ(indn)+4));
                    for i = 1:power
                        subQ{n} = strcat(subQ{n},termStr{n}([indQ(indn)+1,indQ(indn)+2]));
                    end
                else
                    subQ{n} = strcat(subQ{n},termStr{n}([indQ(indn)+1,indQ(indn)+2]));
                end
            end
            
            isZero(n) = rem(count(subQ{n},'1'),2)~=0 || ...
                rem(count(subQ{n},'2'),2)~=0 || rem(count(subQ{n},'3'),2)~=0;
        end
        
        terms = terms(~isZero);
        termStr = termStr(~isZero);
        subQ = subQ(~isZero);
        
        %% combine repeated terms
        % second order
        coeff2 = cell(6+6,1);
        coeff2(:) = {0};
        
        terms2 = terms(cellfun(@(x)length(x)==4,subQ));
        N = length(terms2);
        for n = 1:N
            subQs = zeros(1,4);
            for i = 1:4
                subQs(i) = str2double(subQ{n}(i));
            end
            
            if subQs(1)>subQs(2)
                subQs = [subQs(2),subQs(1),subQs(4),subQs(3)];
            end
            
            if subQs(1)>subQs(3)
                subQs = [subQs(3),subQs(4),subQs(1),subQs(2)];
            end
            
            subQ{n} = num2str(subQs);
            subQ{n} = subQ{n}(~isspace(subQ{n}));
            coeff2ind = cellfun(@(x)strcmp(subQ{n},x),unique2);
            
            indQ = strfind(termStr{n},'Q');
            indN = length(indQ);
            for indn = indN:-1:1
                if strcmp(termStr{n}(indQ(indn)+3),'^')
                    termStr{n}(indQ(indn):indQ(indn)+5) = [];
                else
                    termStr{n}(indQ(indn):indQ(indn)+3) = []; 
                end
            end
            
            coeff2{coeff2ind} = coeff2{coeff2ind}+str2sym(termStr{n});
        end
        
        % third order
        coeff3 = cell(10+18+4,1);
        coeff3(:) = {0};
        
        terms3 = terms(cellfun(@(x)length(x)==6,subQ));
        N = length(terms3);
        for n = 1:N
            subQs = zeros(1,6);
            for i = 1:6
                subQs(i) = str2double(subQ{n}(i));
            end
            
            subii = [subQs(1)==subQs(2),subQs(3)==subQs(4),subQs(5)==subQs(6)];
            if sum(subii)==3
                % iijjkk
                subQs = sort(subQs);
            elseif sum(subii)==1
                % iijkjk
                switch find(subii)
                    case 1
                    case 2
                        subQs = [subQs(3),subQs(4),subQs(setdiff(1:6,[3,4]))];
                    case 3
                        subQs = [subQs(5),subQs(6),subQs(setdiff(1:6,[5,6]))];
                end
                
                if subQs(3)>subQs(4)
                    subQs(3:6) = [subQs(4),subQs(3),subQs(6),subQs(5)];
                end
                
                if subQs(3)>subQs(5)
                    subQs(3:6) = [subQs(5),subQs(6),subQs(3),subQs(4)];
                end
            else
                ascend = [subQs(1)>subQs(2),subQs(3)>subQs(4),subQs(5)>subQs(6)];
                if sum(ascend)<=1
                    subQs = [subQs(2),subQs(1),subQs(4),subQs(3),subQs(6),subQs(5)];
                end
                
                if sum(unique(subQs([1,3,5])))==6
                    subQs = [1,2,2,3,3,1];
                elseif setdiff([1,2,3],[setdiff([1,2,3],unique(subQs([1,3,5]))),...
                        setdiff([1,2,3],unique(subQs([2,4,6])))])==1
                    subQs = [2,1,1,3,2,3];
                elseif setdiff([1,2,3],[setdiff([1,2,3],unique(subQs([1,3,5]))),...
                        setdiff([1,2,3],unique(subQs([2,4,6])))])==2
                    subQs = [1,2,2,3,1,3];
                elseif setdiff([1,2,3],[setdiff([1,2,3],unique(subQs([1,3,5]))),...
                        setdiff([1,2,3],unique(subQs([2,4,6])))])==3
                    subQs = [1,3,3,2,1,2];
                else
                    error('Something unexpected happened');
                end
            end
            
            subQ{n} = num2str(subQs);
            subQ{n} = subQ{n}(~isspace(subQ{n}));
            coeff2ind = cellfun(@(x)strcmp(subQ{n},x),unique3);
            
            indQ = strfind(termStr{n},'Q');
            indN = length(indQ);
            for indn = indN:-1:1
                if strcmp(termStr{n}(indQ(indn)+3),'^')
                    termStr{n}(indQ(indn):indQ(indn)+5) = [];
                else
                    termStr{n}(indQ(indn):indQ(indn)+3) = []; 
                end
            end
            
            coeff3{coeff2ind} = coeff3{coeff2ind}+str2sym(termStr{n});
        end
        
        %% recover symbolic
        if ~isempty(terms2)
            for i = 1:12
                temp = factor(coeff2{i});
                output(n1,n2) = output(n1,n2)+prod(temp)*unique2Sym{i};
            end
        end
        if ~isempty(terms3)
            for i = 1:32
                temp = factor(coeff3{i});
                output(n1,n2) = output(n1,n2)+prod(temp)*unique3Sym{i};
            end
        end
        
        %% generate MATLAB code
        if ~isempty(terms2)
            for i = 1:12
                if coeff2{i}==0
                    continue;
                else
                    coeffStr = coeffMatlab(coeff2{i});
                    MatlabStr{n1,n2} = strcat(MatlabStr{n1,n2},'(',coeffStr,')*',unique2Matlab{i},' + ');
                end
            end
        end
        
        if ~isempty(terms3)
            for i = 1:32
                if coeff3{i}==0
                    continue;
                else
                    coeffStr = coeffMatlab(coeff3{i});
                    MatlabStr{n1,n2} = strcat(MatlabStr{n1,n2},'(',coeffStr,')*',unique3Matlab{i},' + ');
                end
            end
        end
    end
end

end


function [str] = coeffMatlab(coeff)

coeff = char(coeff);
coeff = coeff(~isspace(coeff));
ind = strfind(coeff,'+');
ind = [ind,strfind(coeff,'-')];
ind = [ind,strfind(coeff,'*')];
ind = sort(ind);
if isempty(ind) || ind(1)~=1
    % if the first character is not '-', add a dummy index
    ind = [0,ind];
end
% add a dummy index at end
ind = [ind,length(coeff)+1];

str = '';
N = length(ind)-1;
for n = 1:N
    subCoeff = coeff(ind(n)+1:ind(n+1)-1);
    letterInd = find(isletter(subCoeff));
    
    if isempty(letterInd)
        if n==1 && ind(1)==1
            str = strcat(str,coeff(1),subCoeff);
        else
            str = strcat(str,subCoeff);
        end
    else
        if n==1 && ind(1)==1
            str = strcat(str,coeff(1),subCoeff(letterInd));
        else
            str = strcat(str,subCoeff(letterInd));
        end

        indWedge = strfind(subCoeff,'^');
        if isempty(indWedge)
            str = strcat(str,'(');
            for i = letterInd(end)+1:length(subCoeff)
                str = strcat(str,subCoeff(i),',');
            end
            str(end) = ')';
        else
            str = strcat(str,'(');
            for i = letterInd(end)+1:indWedge-1
                str = strcat(str,subCoeff(i),',');
            end
            str(end) = ')';
            str = strcat(str,subCoeff(indWedge:end));
        end
    end
    
    if n<N
        str = strcat(str,coeff(ind(n+1)));
    end
end

end

