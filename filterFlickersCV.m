function bad_points = filterFlickersCV(data, varargin)
%filters non-gaussian flickers from a data set.
% INPUTS
% data          a 1xN array of time series data points
%
% OPTIONS
% PFLICKER      a 1x2 array of the t-value thresholds for calling a flicker
%               in each of the two stages. default is [3 3].
% WINDOWS       a 1x2 array of the window sizes for finding flickers in
%               each of the two stages. default is [2 5].
%
%
% OUTPUTS
% bad_points    a 1xN logical array which is true where flickers were found


pfs = [3 3];
windows = [2 5];

for ca = 1:length(varargin)
    switch upper(char(varargin{ca}))
        case 'PFLICKER'
            pfs = varargin{ca+1};
        case 'WINDOWS'
            windows = varargin{ca+1};
    end
end

%set up for first stage
bad_data = true;
w = windows(1);
pf = pfs(1);
filtered_data = data;
ii = 1:length(data); %array keeping track of raw indices of filtered data

%loop until nothing is removed
while any(bad_data)
    
    %calculate the median and standard deviation of the points surrounding
    %each data point
    ix = repmat([-w:-1 1:w], length(filtered_data)-2*w, 1) ...
            + repmat((w+1:length(filtered_data)-w)', 1, 2*w);
    data_med = median(filtered_data(ix),2);
    data_std = std(filtered_data(ix),1,2);
    
    %calculate the t-values of the data
    data_deviation = (filtered_data(w+1:end-w)'-data_med)./data_std;
    
    %find data that is outside the threshold
    bad_data = [false(w,1); abs(data_deviation) > pf; false(w, 1)];
    
    %delete that data
    filtered_data(bad_data) = [];
    ii(bad_data) = [];
    
end

%repeat the process with the stage-2 window and t-threshold
bad_data = true;
w = windows(2);
pf = pfs(2);

while any(bad_data)
    ix = repmat([-w:-1 1:w], length(filtered_data)-2*w, 1) ....
            + repmat((w+1:length(filtered_data)-w)', 1, 2*w);
    data_med = median(filtered_data(ix),2);
    data_std = std(filtered_data(ix),1,2);
    data_deviation = (filtered_data(w+1:end-w)'-data_med)./data_std;
    bad_data = [false(w,1); abs(data_deviation) > pf; false(w, 1)];
    filtered_data(bad_data) = [];
    ii(bad_data) = [];
end


%mark removed points as "bad points"
bad_points = true(1, length(data));
bad_points(ii) = false;


end