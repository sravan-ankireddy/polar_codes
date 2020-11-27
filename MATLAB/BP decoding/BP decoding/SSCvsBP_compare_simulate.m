% Length of codeword
N = 256;

% Rate of code
rate = 0.5;

% Length of message vector
K = N*rate;

% Loading Channel indices in decreasing order of reliability
load("indices_" + N + ".mat");
A = true(1,N);
A(ind(K+1:N)) = 0;
A = bitrevorder(A);

% Dense Parity check matrix calculation
H_dense = createDensePolarH(N,A);

% Biparted Parity check matrix calculation
H_sparse = polar2bipartite(N,A,1);
H_sparse = pruneGraph(H_sparse,N);

% Number of message vectors transmitted
num_Tx = 10e2;

% Energy per bit
Eb = 1;

% Eb/N0 step size
delta = 0.3;

% Eb/N0 range of signal in dB
EbN0dB = 1:delta:4;

% Varinace of noise in linear scale
var_N = (10.^(-EbN0dB/10))/(2 * rate);

% Vector to store Bit Error Rate
BER_SSC = zeros(1,length(var_N));
BER_BP = zeros(1,length(var_N));

% Indices corresponding to data channels
data_indices = ind(1:rate*N);
data_indices_sorted = sort(data_indices);

% Indices corresponding to frozen channels
frozen_indices = ind(rate*N + 1:end);
frozen_indices_sorted = sort(frozen_indices);

% Node positions
global Rate_0_node_pos_array;
global Rate_1_node_pos_array;
global SPC_node_pos_array;
global Rep_node_pos_array;
global Rep_SPC_node_pos_array;

Rate_0_node_pos_array = zeros(N);
Rate_1_node_pos_array = zeros(N);
SPC_node_pos_array = zeros(N);
Rep_node_pos_array = zeros(N);
Rep_SPC_node_pos_array = zeros(N);

% Generate node positions
node_position(N, 1, data_indices, data_indices_sorted, frozen_indices, frozen_indices_sorted);

tic;

% Simulation over the range of EbN0
for i_E = 1:length(EbN0dB)

	% Display the simulation status
	disp("counter " + i_E + " of " + length(EbN0dB));
    disp(' ');

	% Std Dev of AWGN noise
	sig = sqrt(var_N(i_E));
    
    %Initialising the number of errors to zero for each sigma
    err_SSC = 0;
    err_BP = 0;

	% Simulation over the message vectors
	parfor i_m = 1:num_Tx

		% Generating a random message vector
        msg = randi([0,1],1,K);

		% Tx vector initialisation
		u0 = zeros(1,N);

		% Assigning data to data channels
		u0(A) = msg;
        
		% Encoding the Tx vector
		xe = encode(u0);

        ub = zeros(N,1);
        ub(A) = msg;
        xb = polarTransform(ub, A);

		% BPSK modulation
		x = sqrt(Eb)*(-1).^xe;

		% Adding AWGN noise
		y = x + sig * randn(size(x));
        yb = xb + sig * randn(size(xb)); 
        
        % LLR Initialisation
		LLR = zeros(log2(N) + 1, N);

		% Channel LLR calculationll
		L_ch = 2*y/(sig^2);
        L_chb = 2*yb/(sig^2);
        LLR(1,:) = L_ch;
 		
        % Simplified Successive Cancellation Decoding
        [LLR, x_hat_SSC, beta] = decode(LLR, 1, 1, Rate_0_node_pos_array, Rate_1_node_pos_array, SPC_node_pos_array, Rep_node_pos_array, Rep_SPC_node_pos_array);
        
        % Belief Propagation decoding
        max_iter = 200;
        Lch_ext=zeros(size(H_sparse,2),1);
        Lch_ext((end-N+1):end)=L_chb;
        [msg_hat, ~, ~] = ldpc_minsum_decode(H_sparse,0,Lch_ext,max_iter);
        msg_hat = msg_hat((end-N+1):end);
        msg_hat = polarTransform(msg_hat, A);
        msg_hat = msg_hat(A);

        % Checking for errors
        err_SSC = err_SSC + sum(xor(u0(data_indices),x_hat_SSC(data_indices)));
        err_BP = err_BP + sum(xor(msg,msg_hat));
                
    end    
    % Bit Error Rate
    BER_SSC(1,i_E) = err_SSC/(num_Tx*N);
    BER_BP(1,i_E) = err_BP/(num_Tx*N);
end

disp(BER_SSC);
disp(BER_BP);

toc;
%% plot for BER_SSCL vs Eb/N0
figure(1)
semilogy(EbN0dB,BER_SSC,'-s');
hold on;
semilogy(EbN0dB,BER_BP,'-p');

legend('SSC', 'BP');
str = sprintf('Plot for BER SSC, BP vs Eb/N0 for length L = %d, %d Tx vectors', N, num_Tx);
title(str);
xlabel('Eb/N0 in dB');
ylabel('Bit Error Rate(BER)');
hold on; 