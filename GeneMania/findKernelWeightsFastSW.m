function [indices,alpha] = findKernelWeightsFastSW(y,W)
% 	function [indices,alpha] = findKernelWeightsFastSW(y,W)

[NN,CC] = size(y);
MM = length(W);
Nij = sparse([],[],[],NN,NN,0);

epsilon = eps;             % our numerical def'n of zero

% --- Calculation of number of time i is positive and number ot times (i,j)
% is positive
[Nij,BHalf,const] = fastcoann(y);

Nij = Nij - diag(diag(Nij)); % subtract counting of (i,i)
Ni = sum(y,2);

biasVal = NN*NN*CC;
KtK = zeros(MM+1,MM+1);
KtT = zeros(MM+1,1);
KtK(1,1) = biasVal;
KtT(1) = sum(BHalf'*NN)+sum(sum(Nij))+const*NN*NN;

for ii = 1:MM
    KtK(ii+1,1) = CC*sum(sum(W{ii}));
    KtK(1,ii+1) = KtK(ii+1,1);
    KtT(ii+1,1) = sum(sum(W{ii}.*Nij))+sum(BHalf'*W{ii})-sum(BHalf'*(W{ii}.*Nij))+sum(sum(W{ii}))*const;
   
    for jj = 1:ii
        KtK(ii+1,jj+1) = CC*sum(sum(W{ii}.*W{jj}));
        KtK(jj+1,ii+1) = KtK(ii+1,jj+1);
    end
end

KtK = full(KtK);
KtT = full(KtT);


% ---- Removing empty columns

ss = sum(abs(KtK), 2);
goodIdx = find(sum(abs(KtK), 2) > epsilon * max(ss));
KtK = KtK(goodIdx, goodIdx);
KtT = KtT(goodIdx);

indices = 1:(length(goodIdx) - 1);

done = 0;
if nargin > 2
    oneI = 1*eye(size(KtK));
    oneI(1,1) = 0;
    KtK = KtK +oneI*100;
end

while ~done,

    % alpha = (K' * K)^-1 K' * T(:)
    alpha = inv(KtK) * KtT;
    if(any(isnan(alpha)))
        indices = [];
        break
    end

    % Find the locations of all the negative weights
    negWeights = find(alpha < 0);

    negWeights = negWeights(:);
    negWeights = setdiff(negWeights, 1);  % always have a bias term now
    if length(negWeights) > 0,

        % Save ourselves some computation by just eliminating rows and cols
        % from K^T K and elements from K^T y. The savings from this appear
        % to be nominal, though.
        KtK(:,negWeights) = [];
        KtK(negWeights,:) = [];
        KtT(negWeights) = [];

        % Remove from the indices array.
        indices((negWeights - 1)) = [];
    else,
        done = 1;
        fprintf('%d matrices chosen, indices %s, weights %s\n',...
            size(KtK,1)-1,mat2str(indices),mat2str(alpha(2:end)));
    end;
end;




% Undo our previous trickery
indices = goodIdx(indices);


% if we have no non-zero alpha entries, use the average kernel
if length(indices) == 0,
  indices = 1:length(W);
  alpha = [0 1/length(W)*ones(1,length(W))]';
  fprintf('All kernels eliminated or empty, using average kernel\n');
end;
% ignore the bias term which is alpha(1)
alpha = alpha(2:end);
        
    
    



