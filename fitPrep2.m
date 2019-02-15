function [x1,x2,ix] = fitPrep2(x1,x2)
%Formats two sets of data which may or may not be column vectors and may
%or may not contain nans into fit-appropriate nan-free column vector

if size(x1,1) < size(x1,2)
    x1 = x1';
end

if size(x2,1) < size(x2,2);
    x2 = x2';
end

ix = find(~isnan(x1) & ~isnan(x2));
x1 = x1(ix); x2 = x2(ix);

end