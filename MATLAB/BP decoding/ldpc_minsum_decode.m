function [rcap, Ltot, it, invalid]=ldpc_minsum_decode(H,offset,L0,maxitrs)

L0 = L0';

%Setup
[row,col]=find(H); [nrow,n]=size(H);
ind=sub2ind([nrow,n],row,col);

%Decoding
Ltot = L0;
L=L0(col); rcap=(L0<0);
it=0; invalid=any(bitand(H*rcap',1));

while ((it<maxitrs)&&(invalid))
    s=sparse(row,col,L<0);
    P=bitand(full(sum(s,2)),1);

    L=full(abs(L));
    for i=1:nrow
        indi=(row==i);
        Li=L(indi);
        [mini,loc]=min(Li);
        min2i=min(Li([1:loc-1 loc+1:end]));
        Li(1:end)=mini-offset;
        Li(loc)=min2i-offset;
        L(indi)=Li;
    end
    
    L=sparse(row,col,(1-2*P(row)).*(1-2*s(ind)).*L');
    
    Ltot=L0+full(sum(L)); rcap=(Ltot<0);
    if any(bitand(H*rcap',1))
        L=Ltot(col)-L(ind');
    else
        invalid=0;
    end
    it=it+1;
    
end


