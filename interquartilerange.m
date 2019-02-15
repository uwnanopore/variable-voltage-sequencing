function [iqr, plus_iqr, minus_iqr, iqr_points] = interquartilerange(x)


if isempty(x)
    iqr = nan;
    return;
end

x = sort(x, 2);

iqr_indices = (0:4)*(size(x,2)-1)/4+1;
iqr_index_fractions = ones(size(x, 1),1)*mod(iqr_indices,1);
iqr_lower_indices = floor(iqr_indices);
iqr_upper_indices = ceil(iqr_indices);

iqr_points = x(:,iqr_lower_indices).*iqr_index_fractions + x(:,iqr_upper_indices).*(1-iqr_index_fractions);
iqr = iqr_points(:,4) - iqr_points(:,2);
plus_iqr = iqr_points(:,4) - iqr_points(:,3);
minus_iqr = iqr_points(:,3) - iqr_points(:,2);


end