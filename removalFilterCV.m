function filteredread = removalFilterCV(read)

bounds= [.10,.85];
mindur = 20;



goods = (read.npts_uf >= mindur & read.x_uf(1,:) >= bounds(1) & read.x_uf(1,:) <= bounds(2));
%update jumps
ix = 1:length(read.x_uf);
ix(~goods) = nan;
read.x_f = read.x_uf(:,goods);
read.k_f = read.k_uf(goods);
read.uf_to_f = ix;

filteredread =read;

end