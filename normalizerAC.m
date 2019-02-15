function padded_n = normalizerAC(features)

good_indices = ~isnan(features(1,:));
padded_n = nan(size(features));
features = features(:,good_indices);
%find the average curve
n0 = nanmean(features, 2);

%find the normalized curves
g = features - n0;
gmeans = mean(g,1);

%indices for voltage points
v = ((1:size(g,1)) - mean(1:size(g,1)))';


m = zeros(1, size(features, 2));
for ii = 1:length(gmeans)
    p = polyfit(v, g(:,ii), 1);
    m(ii) = p(1);
end

[gmeans, m] = fitprep(gmeans, m);
pm = polyfit(gmeans, m, 1);


n = zeros(size(features));
for ii = 1:length(gmeans)
    n(:,ii) = n0 + (pm(1)*g(:,ii) + pm(2)).*v;
end

padded_n(:,good_indices) = n;

end