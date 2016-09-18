function net = load_network_text( file_name, net_g2i)
%LOAD_NETWORK_TEXT Summary of this function goes here
%   Detailed explanation goes here

   [e1,e2,wt,~] = textread(file_name,'%s%s%f%s');
      e1 = upper(e1);
   e2 = upper(e2);
   e1 = cell2mat(values(net_g2i,e1));
   e2 = cell2mat(values(net_g2i,e2));
   [~,IA,~] = unique([e1,e2],'rows');
   new_e1 = e1(IA);
   new_e2 = e2(IA);
   new_wt = wt(IA);
   ngene = size(net_g2i,1);
   net = sparse(new_e1,new_e2,new_wt,ngene,ngene);   
   assert(max(net(:))<=1);
end

