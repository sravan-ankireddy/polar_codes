% close all;
clear all;

% Length of codeword
N = 1024;

% Loading Channel indices in decreasing order of reliability
load("indices_" + N + ".mat");

% Rate of code
rate = 0.5;

% Length of message vector
K = N*rate;

% Number of message vectors transmitted
num_Tx = 10e3;

% Energy per bit
Eb = 1;

% Eb/N0 step size
delta = 0.3;

% Eb/N0 range of signal in dB
EbN0dB = 1:delta:4;
% EbN0dB = 1;

% Varinace of noise in linear scale
var_N = (10.^(-EbN0dB/10))/(2 * rate);

% Vector to store Bit Error Rate
BER_SSC = zeros(1,length(var_N));

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
    err_SSCL = 0;
    err_SCL = 0;
    err_SSC = 0;

	% Simulation over the message vectors
	parfor i_m = 1:num_Tx

		% Generating a random message vector
        msg = randi([0,1],1,K);

		% Tx vector initialisation
		u0 = zeros(1,N);

		% Assigning data to data channels
		u0(data_indices) = msg;
        
		% Encoding the Tx vector
		xe = encode(u0);

		% BPSK modulation
		x = sqrt(Eb)*(-1).^xe;

		% Adding AWGN noise
		y = x + sig * randn(size(x));

		% SSCL Decoding at the reciever
        in = 1;
        ind_lev = 1;
        
        % LLR Initialisation
		LLR = zeros(log2(N) + 1, N);

		% Channel LLR calculationll
		LLR(1,:) = 2*y/(sig^2);
 		
        [LLR, y_hat_SSC, beta] = decode(LLR, in, ind_lev, Rate_0_node_pos_array, Rate_1_node_pos_array, SPC_node_pos_array, Rep_node_pos_array, Rep_SPC_node_pos_array);
        
        % Checking for errors
        err_SSC = err_SSC + sum(xor(u0(data_indices),y_hat_SSC(data_indices)));
                
    end    
    % Bit Error Rate
    BER_SSC(1,i_E) = err_SSC/(num_Tx*N);
    
end

disp(BER_SSC);

toc;
%% plot for BER_SSCL vs Eb/N0
figure(1)
semilogy(EbN0dB,BER_SSC,'-s');
hold on;

legend('SSC', 'SCL', 'SSCL');
% str = sprintf('Plot for BER SSC, SCL, SSCL vs Eb/N0 for length %d and L = %d, %d Tx vectors', N, l, num_Tx);
% title(str);
xlabel('Eb/N0 in dB');
ylabel('Bit Error Rate(BER)');
hold on; 