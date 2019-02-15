function varargout = fitprep(varargin)
%Formats sets of data which may or may not be column vectors and may or
%may not contain nans into column vectors with nans removed for use with the
%function fit().
%
%inputs can either be vectors for formatting, or one cell array of vectors
%for formatting.
%
%outputs are, correspondingly, either vectors that have been formatted or a
%cell array of formatted vectors. the last output is always the indices of
%the data that were passed through.

if iscell(varargin{1})
    data = varargin{1};
else
    data = varargin;
end

for cd = 1:length(data)
    if size(data{cd},1) < size(data{cd},2)
        data{cd} = data{cd}';
    end
    if isempty(data{cd})
        varargout = cell(1, length(data));
        return
    end
end
data = cell2mat(data);

ix = find(~any(isnan(data) | isinf(data), 2));
data = data(ix,:);
if iscell(varargin{1})
    varargout{1} = num2cell(data, 1);
else
    varargout = num2cell(data, 1);
end
varargout{end+1} = ix;

end