function [X, means, stddevs] = whiten(X, means, stddevs)

if ~exist('means');
    means = nanmean(X);
end

if ~exist('stddevs')
    stddevs = nanstd(X);
end


[n p] = size(X);

X = bsxfun(@minus, X, means);

X(find(isnan(X))) = 0;

X = bsxfun(@rdivide, X, stddevs);

X(find(isnan(X))) = 0;
