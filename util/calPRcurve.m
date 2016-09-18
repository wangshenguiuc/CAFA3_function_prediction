function [rc,pr,r,u]= calPRcurve(pred,ref,tau,w)

  k = numel(tau);

  nT = sum(ref); % the number of positives in the reference
  
  pos = ref;
  neg = ~ref;
  

  wpos = sum(w .* pos);
  wneg =  sum(w .* neg);
  
  pr = zeros(k, 1);
  rc = zeros(k, 1);
  
  r = zeros(k, 1);
  u = zeros(k, 1);
  for i = 1: k
    P  = (pred >= tau(i));
    nP = sum(P); % the number of each predicted positives

    nTP =sum(P & ref);
    pr(i) = nTP ./ nP;
    rc(i) = nTP ./ nT;
    
    r(i) = dot(w,(P & neg));
    u(i) =  wpos - dot(w,(P & pos));
  end

  
end
