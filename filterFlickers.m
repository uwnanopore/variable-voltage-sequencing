function [bad_cycles] = filterFlickersAC(data, varargin)
% filter non-gaussian behaviour from variable-voltage data

pfs = [3, 3]; % the score thresholds for point removal
windows = [2, 5]; % the window sizes in number of points

for ca = 1:length(varargin)
    switch upper(char(varargin{ca}))
        case 'PFLICKER'
            pfs = varargin{ca+1};
        case 'WINDOWS'
            windows = varargin{ca+1};
    end
end

%a list of bad points. initialize to true to allow loop to begin
bad_data = true;

%use the stage-one parameters
w = windows(1);
pf = pfs(1);

%the filtered product
filtered_data = data;

%indexing to track where the preserved points were in the original data-
%this is what is returned by the whole function
ii = 1:length(data);

%repeat until no changes are made
while any(bad_data)
    
    %indices of the points in each window surrounding each point, arranged
    %in an N x 2w array; each row is one window corresponding to the point
    %at its center.
    ix = repmat([-w:-1 1:w], length(filtered_data)-2*w, 1) ...
                    + repmat((w+1:length(filtered_data)-w)', 1, 2*w);
                
    %find the median and standard deviation of the data in each window
    data_med = median(filtered_data(ix),2);
    data_std = std(filtered_data(ix),1,2);
    
    %compute the t-value for the central point
    data_deviation = (filtered_data(w+1:end-w)'-data_med)./data_std;
    
    %eliminate data whose t value is sufficiently large to indiciate it is
    %an outlier from the surrounding points
    bad_data = [false(w,1); abs(data_deviation) > pf; false(w, 1)];
    filtered_data(bad_data) = [];
    ii(bad_data) = [];
    
end


%repeat the above process but with more aggressive parameters
bad_data = true;
w = windows(2);
pf = pfs(2);

while any(bad_data)
    ix = repmat([-w:-1 1:w], length(filtered_data)-2*w, 1) ...
                    + repmat((w+1:length(filtered_data)-w)', 1, 2*w);
    data_med = median(filtered_data(ix),2);
    data_std = std(filtered_data(ix),1,2);
    data_deviation = (filtered_data(w+1:end-w)'-data_med)./data_std;
    bad_data = [false(w,1); abs(data_deviation) > pf; false(w, 1)];
    filtered_data(bad_data) = [];
    ii(bad_data) = [];
end


bad_cycles = true(1, length(data));
bad_cycles(ii) = false;


end