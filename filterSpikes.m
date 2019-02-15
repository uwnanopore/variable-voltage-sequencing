function [x_filtered, filtered_indices] = filterSpikes(x, varargin)

score_cutoff = 5;
max_removed = 10;
for ca = 1:length(varargin)
    switch upper(char(varargin))
        case 'CUTOFF'
            score_cutoff = varargin{ca+1};
        case 'MAXREMOVED'
            max_removed = varargin{ca+1};
    end
end

% load AC_princomps_normalized
pc = load('principal_components.mat');
pc = pc.principal_components;
x_filtered = x;
filtered_indices = false(size(x));
for cl = 1:size(x,2)
    
    goods = true(size(x,1),1);
    num_filtered = inf;
    while num_filtered > 0
        A = pc(goods,:);
        smoother = A*((A'*A)\A');
        this_x = x(:,cl);
        ix = find(goods);
        resids = this_x(goods) - smoother*this_x(goods);
        scores = resids.^2/var(resids);
        bads = scores > score_cutoff;
        num_filtered = sum(bads);
        goods(ix(bads)) = false;
    end
    
    A = pc(goods, :);
    B = pc(~goods, :);
    if sum(~goods) <= max_removed
        x_filtered(~goods,cl) = B*((A'*A)\A')*this_x(goods);
        filtered_indices(:,cl) = ~goods;
    end
end



end