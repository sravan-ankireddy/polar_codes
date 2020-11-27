function [rcap, L1, it, invalid]=ldpc_minsum_decode_alternate(H_BP_tot,H_par,N,num_par_eq,offset,L0,maxitrs)

L0 = L0';

L10 = zeros(size(H_BP_tot,2),1); % belongs to H_BP
L20 = zeros(size(H_par,2),1); % belongs to H_par


L10 = L10';
L20 = L20';

%Setup
[row_BP,col_BP]=find(H_BP_tot); [nrow_BP,n_BP]=size(H_BP_tot);
ind_BP=sub2ind([nrow_BP,n_BP],row_BP,col_BP);

[row,col]=find(H_par); [nrow,n]=size(H_par);
ind=sub2ind([nrow,n],row,col);

%Decoding
L10(end-N+1:end) = L0(1:N);
L1=L10(col_BP); rcap=(L10<0);

L20 = L0;
L0a = L0;
L2=L20(col);
% whos H_BP; whos rcap; whos L1;
it=0; invalid=any(bitand(H_BP_tot*rcap',1));
while ((it<maxitrs)&&(invalid))
    whos row_BP; whos col_BP; whos L1;
    
    s_BP=sparse(row_BP,col_BP,L1<0);
    P_BP=bitand(full(sum(s_BP,2)),1);
    
    Ls = L1;
    Ltot_old=L10+full(sum(L1));
    L1=full(abs(L1));
    for i=1:nrow_BP
        indi=(row_BP==i);
        Li=L1(indi);
        [mini,loc]=min(Li);
        min2i=min(Li([1:loc-1 loc+1:end]));
        Li(1:end)=mini-offset;
        Li(loc)=min2i-offset;
        L1(indi)=Li;
    end

    L1=sparse(row_BP,col_BP,(1-2*P_BP(row_BP)).*(1-2*s_BP(ind_BP)).*L1');

    Ltot=L10+full(sum(L1)); rcap=(Ltot<0);
%     whos Ltot; whos L10; whos rcap;
    if any(bitand(H_BP_tot*rcap',1))
        L1=Ltot(col_BP)-L1(ind_BP');
    else
        invalid=0;
    end
    whos L2;
    if (it == 0)
        L0a(1:N) = Ltot(end-N+1:end);
    else
        L0a(1:N) = L0(1:N) + Ltot(end-N+1:end) - Ltot_old(end-N+1:end);% - Ls(end-N+1:end);
    end
    whos H_par; whos row; whos col; whos L2a;
    L2=L0a(col);
    s=sparse(row,col,L2<0);
    P=bitand(full(sum(s,2)),1);
        
    Le = L2;
    Ltot_old = L20+full(sum(L2));
    L2=full(abs(L2));
    for i=1:nrow
        indi=(row==i);
        Li=L2(indi);
        [mini,loc]=min(Li);
        min2i=min(Li([1:loc-1 loc+1:end]));
        Li(1:end)=mini-offset;
        Li(loc)=min2i-offset;
        L2(indi)=Li;
    end
    whos L2;
    L2=sparse(row,col,(1-2*P(row)).*(1-2*s(ind)).*L2');
    whos L2;
    Ltot=L0+full(sum(L2)); %rcap=(Ltot<0);
%     if any(bitand(H_BP*rcap',1))
        L2=Ltot(col)-L2(ind');
%     else
%         invalid=0;
%         break;
%     end
%     whos L1; whos L0; whos L2; whos Le;
    L0a(1:N) = L0(1:N) + Ltot(1:N) - Ltot_old(1:N);% - Le(1:N);
    it=it+1;
    
end


