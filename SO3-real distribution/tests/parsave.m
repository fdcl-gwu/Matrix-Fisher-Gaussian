function [  ] = parsave( filename, varargin )

for i = 1:length(varargin)
    temp.(inputname(i+1)) = varargin{i};
end

save(filename,'-struct','temp','-v7.3');

end

