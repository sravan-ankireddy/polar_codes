EbNodB_mat = 1:0.3:4;
 Nblocks = 100;
BER = zeros(length(EbNodB_mat),1);
FER = zeros(length(EbNodB_mat),1);

tic;
for i = 1:length(EbNodB_mat)
    disp(i);
    EbNodB = EbNodB_mat(i);
    nrpolar_sclistdecode_FP;
    BER(i) = BER_sim;
    FER(i) = FER_sim;
end
toc;

figure(1)
semilogy(EbNodB_mat, BER);
hold on;
semilogy(EbNodB_mat, FER);
