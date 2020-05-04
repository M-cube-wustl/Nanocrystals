function [ out ] = dircat( varargin )
%concatonate character arrays with '\' in between 

out = varargin{1};
N = numel(varargin);

for i=2:N
    out = strcat(out,'\',varargin{i});
end
 
end

