function ye = encode(ue)
N = length(ue);
n = log2(N);

ye = ue;

for i_s = 1:n                       %stage index
    for i_g = 1:2^(n-i_s)           %group index
        del = 2^(i_s-1);            %delta
        base = (2^i_s)*(i_g - 1);   %base for sub-group
        for i_sg = 1:2^(i_s-1)      %sub-group index
          ye(base+i_sg) = xor(ye(base+i_sg), ye(base+i_sg+del));
        end
         
    end
end
end
