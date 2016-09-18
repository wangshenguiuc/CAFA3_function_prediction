function [coAnn,BHalf,const] = fastcoann(Y);

[N,D] = size(Y);
N2 = N*N;
% calculate category sizes
sizes = full(sum(Y, 1));  %% have to make sure that all values are 0 or 1
cache = cell(1,max(sizes));       %% caches the two-column index matrices made by allpairs.m
cachePN = cell(1,max(sizes));



%%  Calculate how much memory we will need to store the elements of the
%%  sparse matrix
veclen = 0; 
veclenPN = 0;
for ii = 1:length(sizes); 
    veclen = veclen + (sizes(ii))*(sizes(ii)-1) / 2; 
    veclenPN = veclenPN + sizes(ii);
end

%% pre-allocate vectors used to make the sparse matrix
rows = zeros(1,veclen);   % will hold the row indices for the co-annotated pair entries
cols = zeros(1,veclen);   % will hold the col indices
vals = zeros(1,veclen);   % will hold the constant

rowsPN = zeros(1,veclenPN);
colsPN = zeros(1,veclenPN);
valsPN = zeros(1,veclenPN);

% Note that we only populate one half of the co-annotation matrix, we'll
% add it to its transpose afterwards

offset = 0;  % stores where we are in the vectors 
offsetPN = 0;
const = 0;
for ii = 1:length(sizes); 
if mod(ii,10000)==0
fprintf('finished %f\n',ii/length(sizes));
end
    if sizes(ii) > 0; 
        nn = sizes(ii); 
        if isempty(cache{nn})  % check if we've already calculated the two column index matrix 
            cache{nn} = allpairs(nn);  % if not, calculate and cache it.
        end 
        if isempty(cachePN{nn});
            cachePN{nn} = 1:nn;
        end
        pp = cache{nn};         % pp is the index matrix
        pn = cachePN{nn};
        len = size(pp, 1);      % # of elements we are added to the vectors
        
        lenPN = length(pn);
        
        myrows = find(Y(:,ii));   % row indices of the annotated genes
        rows(offset+1:offset+len) = myrows(pp(:,1));   % copy row indices
        cols(offset+1:offset+len) = myrows(pp(:,2));   % copy col indices
        vals(offset+1:offset+len) = 1;   % replace this with correct constant 
        
        rowsPN(offsetPN+1:offsetPN+lenPN) = myrows(pn);
        valsPN(offsetPN+1:offsetPN+lenPN) = (-2*nn)./N;
        
        offset = offset+len; 
        offsetPN = offsetPN+lenPN;
        
        const = const + (nn*nn)./(N2);
        
        if mod(ii, 100) == 0;   % status update for inpatient users, remove if not inpatient
            fprintf('.'); 
        end
    end
end

ngenes = size(Y, 1);
% Note sparse adds together elements with the same row and column index, below we rely on this behaviour,
% if re-implementing make sure to do this addition.
lowerHalf = sparse(rows, cols, vals, ngenes, ngenes);  
BHalf = sparse(rowsPN,1,valsPN,ngenes,1);
% As mentioned above, we only populate one of the (i,j) elements for each
% class, here we ensure that the (j,i) element has the same value.  It's
% important to add because we may populate (i,j) for some classes and (j,i)
% for others (though probably not).
coAnn = lowerHalf + lowerHalf';




