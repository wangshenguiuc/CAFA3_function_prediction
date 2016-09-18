function [r,u]= calRUcurve(pred,ref,tau)

  k = numel(tau);

  pos = ref;
  neg = ~ref;
  
  npos = sum(pos);
  nneg = sum(neg);
  
  r = zeros(k, 1);
  u = zeros(k, 1);
  
  for i = 1 : k
    P  = (pred >= tau(i));
    r(i) =sum(P & neg);
    u(i) = npos -  sum(P & pos);
  end
  
  
  
end
