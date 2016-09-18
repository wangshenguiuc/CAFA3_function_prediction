function p = allpairs(n)
%function p = allpairs(n)
%
%  N is a whole number
%  P is a two-column matrix of all sets of size 2 from the set {1, 2, ..., N}

p = zeros(n * (n - 1) / 2, 2);
row = 1;
for ii = 1:n
    for jj = 1:(ii-1)
        p(row, :) = [ii jj];
        row = row + 1;
    end
end
