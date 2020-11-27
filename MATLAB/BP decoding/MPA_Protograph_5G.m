% Rate
rate = 0.5;

% Codeword blocklength
N = 1944;
K = N*rate;

% Subblock length
Z = 81;

% Parity Check Matrix
H = zeros(N - K,N);

% Select base matrix from 5G specs based on rate
H_base = base_matrix(rate);

disp("Calculating Parity Check matrix ... ")
% Expanding the base matrix to obtain final parity check matrix
for i_P = 1:size(H_base,1)
    for j_P = 1:size(H_base,2)
        if (H_base(i_P, j_P) >= 0)
            H( (i_P - 1)*Z + 1: i_P*Z, (j_P - 1)*Z + 1: j_P*Z ) = cyclic_shift(H_base(i_P, j_P), eye(Z));
        end
    end
end
disp(" .. finished");

% finding generator matrix by reducing H to row echelon form
% R = rref(H);
% P = transpose(R(:,N-K+1:end));
% G = [eye(K) P];

% offset for minsum
offset = 0;

% No. of samples in EBN0 range
nSam = 10;

% start and stop of EbN0dB
EbN0dB_start = 1; EbN0dB_stop = 3.7;

% Eb/N0 step size
delta = (EbN0dB_stop - EbN0dB_start)/(nSam-1);

% Eb/N0 range of signal in dB
EbN0dB = EbN0dB_start:delta:EbN0dB_stop;

% Varinace of noise in linear scale
var_N = (10.^(-EbN0dB/10))/(2 * rate);

% no. of simulations
Nblocks = 100*1000;

% no. of iterations
max_iter = 20;

[ParRI, ParCI] = find(H == 1);
ParI = [ParRI, ParCI];

decoder = comm.LDPCDecoder('ParityCheckMatrix',ParI,'OutputValue','Whole codeword','DecisionMethod','Soft decision','MaximumIterationCount',max_iter,'NumIterationsOutputPort',1,'IterationTerminationCondition', 'Parity check satisfied');

% error variables
BER_original = zeros(size(EbN0dB));
BER_LDPC = zeros(size(EbN0dB));
BLER_LDPC = zeros(size(EbN0dB));

tic;
for i_EbN0dB = 1:length(EbN0dB)
    
    disp(i_EbN0dB);
    err_original = 0;
    err_LDPC = 0;
    bler_LDPC = 0;
    err_and = 0;
    parfor i_sim = 1:Nblocks

        % Generating a random message vector
        msg = randi([0,1],K,1);
        
        
        % generate all zero codeword
        c = zeros(N,1);
        
        % 
        
        % BPSK modulation
        x = 1 - 2*c;

        % Adding AWGN noise
        sig = sqrt(var_N(i_EbN0dB));
        y = x + sig * randn(size(x));

        % Rx LLRs
        L_ch = 2*y/(sig^2);
        
        % est before LDPCoding
        c_hat_original = ones(N,1);
        c_hat_original(L_ch > 0) = 0; 

        % bit error rate
        err_original = err_original + sum(c_hat_original ~= c);

        L = step(decoder, L_ch); c_hat_LDPC=(-sign(L))/2+0.5;

        % bit error rate
        err_LDPC = err_LDPC + sum(c_hat_LDPC ~= c);
        if (sum(c_hat_LDPC ~= c) > 0)
            bler_LDPC = bler_LDPC + 1;
        end

    end
    BER_original(i_EbN0dB) = err_original/(N*Nblocks);
    BER_LDPC(i_EbN0dB) = err_LDPC/(N*Nblocks);
    BLER_LDPC(i_EbN0dB) = bler_LDPC/(Nblocks);
end
toc;

disp(BER_original);
disp(BER_LDPC);

%% Plots
figure(1)
% semilogy(EbN0dB, BER_original);
semilogy(EbN0dB, BER_LDPC, '--*');
hold on;
semilogy(EbN0dB, BLER_LDPC, '-*');
xlabel({'Eb/N0', '(in dB)'});
ylabel('BER/BLER');
legend('BER', 'BLER');
title('BER vs Eb/N0');