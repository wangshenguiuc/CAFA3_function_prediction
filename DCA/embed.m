function w= embed(train_ind,A,LX,LY,weight)
    nlabel = size(LY,1);
    nnode = size(LX,1);
    xy=zeros(nnode,nlabel);
    for i=1:nlabel
%         p = intersect(pos_ind(i,1:pos_l(i)),train_ind);
%         n = intersect(neg_ind(i,1:neg_l(i)),train_ind);
%         p1 = intersect(find(A(:,i)==1),train_ind);
%         n1 = intersect(find(A(:,i)==0),train_ind);
%         p = train_ind(A(train_ind,i)==1);
%         n = train_ind(A(train_ind,i)==0);
         p = train_ind(A(train_ind,i)==1);
        n = train_ind(A(train_ind,i)==0);
%         assert(isequal(sort(p1),sort(p'))&&isequal(sort(n),sort(n1')));
%         xy(p,i)=length(n);
%         xy(n,i)=-length(p);
    xy(p,i)=sum(weight(n));
        xy(n,i)=-sum(weight(p));
        
    end
%     spec_nnode = size(specie_LX,1);
    w = LX'*xy*LY;
end