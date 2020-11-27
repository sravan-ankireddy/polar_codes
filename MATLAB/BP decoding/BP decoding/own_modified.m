
tic;
%Simulation Parameters

N = 64;
R = 0.5;
desSNR = 0.6;
max_iter = 200;
nbCWs = 10e3;
nbSNR = 15;
iter = zeros(nbSNR,1);
EBN0dB_start = 1;
EbN0dB_stop = 6;
bipartite = 0;        %use LDPC-like code construction as presented in the paper (naive dense matrix is used instead)

%Initialize variables and
EbN0dB = linspace(EBN0dB_start,EbN0dB_stop,nbSNR);
K = round(N*R);

BER_BP = zeros(nbSNR,1);
BER_SSC = zeros(nbSNR,1);

BLER_BP = zeros(nbSNR,1);
BLER_SSC = zeros(nbSNR,1);

load("indices_" + N + ".mat");
A = true(1,N);
A(ind(K+1:N)) = 0;
A = bitrevorder(A);

% Indices corresponding to data channels
data_indices = ind(1:R*N);
data_indices_sorted = sort(data_indices);

% Indices corresponding to frozen channels
frozen_indices = ind(R*N + 1:end);
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

%%%%%%%%%start pruning%%%%%%%%%%%%%%5
disp('Create and prune H');
%create H matrix

if bipartite==1
    [H] = polar2bipartite(N,A,1); %already kicks out frozen indices
    %prune all useless nodes
    H = pruneGraph(H,N);
else
    H = createDensePolarH(N,A);
    Hdense = H; 
end


%%%%%%%%%%%%START BER SIMULATIONS%%%%%%%%%%%%%%%%%
%initialize decoder
disp('Starting BER simulation');

[ParRI, ParCI] = find(H == 1);
ParI = [ParRI, ParCI];

decoder = comm.LDPCDecoder('ParityCheckMatrix',ParI,'OutputValue','Whole codeword','DecisionMethod','Soft decision','MaximumIterationCount',max_iter,'NumIterationsOutputPort',1,'IterationTerminationCondition', 'Parity check satisfied');

for j=1:nbSNR
    disp(['evaluating SNRb=' num2str(EbN0dB(j)) 'dB, ' num2str(j) '/' num2str(nbSNR)]);
    
    sigma = sqrt((10.^(-EbN0dB(j)/10))/(2 * R));
    %initialize temp variables
    BER_BP_temp = zeros(nbCWs,1);
    BER_SSC_temp = zeros(nbCWs,1);
    nbBit_temp = zeros(nbCWs,1);
    
    BLER_BP_temp = zeros(nbCWs,1);
    BLER_SSC_temp = zeros(nbCWs,1);
    nbBlocks_temp = zeros(nbCWs,1);
    numiter = zeros(nbCWs,1);
    
    parfor i=1:nbCWs
        
        % Generating a random message vector
        msg = randi([0,1],K,1);   
        
        % Tx vector initialisation
		u_BP = zeros(N,1);
        
        u_BP(A) = msg;
        
        x_BP = encode(u_BP);
        
        % BPSK modulation
        x_BP = (1 - 2*x_BP);

		% Assigning data to data channels
        u_SSC = zeros(N,1);
        
		u_SSC(data_indices) = msg;
        
        x_SSC = encode(u_SSC);
        
        x_SSC = (1 - 2*x_SSC);
         
        % Adding AWGN noise
        noise = sigma*randn(size(x_BP));
        
        y_BP = x_BP + noise;
        
        y_SSC = x_SSC + noise;
        
        % Channel LLRs
        Lch = 2*y_BP./sigma.^2;
        Lch_SSC = 2*y_SSC./sigma.^2; 
        
        % extend Lch vector by 0 positions
        Lch_ext = zeros(size(H,2),1);
        xhat_BP = zeros(size(H,2),1);
        
        % assuming the last positions are channel positions
        Lch_ext((end-N+1):end) = Lch;     
                
        [LLRxhat, iter] = step(decoder, Lch_ext); xhat_BP=(-sign(LLRxhat))/2+0.5;    %and decode

        LLR = zeros(log2(N) + 1, N);
        LLR(1,:) = Lch_SSC;
        [~, x_hat_SSC, ~] = decode(LLR, 1, 1, Rate_0_node_pos_array, Rate_1_node_pos_array, SPC_node_pos_array, Rep_node_pos_array, Rep_SPC_node_pos_array);        
        
        msg_hat_BP = xhat_BP((end-N+1):end);
        
        msg_hat_BP = polarTransform(msg_hat_BP, A);
        
        msg_hat_BP = msg_hat_BP(A);
        
        msg_hat_SSC = x_hat_SSC(data_indices)';
                
        %count errors
        BER_BP_temp(i) = sum(msg_hat_BP ~= msg);
        nbBit_temp(i) = length(msg_hat_BP);
        
        BER_SSC_temp(i) = sum(msg_hat_SSC ~= msg);
        
        %also consider BLER
        if sum(msg_hat_BP ~= msg)~=0
            BLER_BP_temp(i)=1;
        end
        
        if sum(msg_hat_SSC ~= msg)~=0
            BLER_SSC_temp(i) = 1;
        end
        nbBlocks_temp(i) = 1;
    end
    %sum up temp arrays (due to parfor loop)
    BER_BP(j) = BER_BP(j)+sum(BER_BP_temp);
    BLER_BP(j) = BLER_BP(j)+sum(BLER_BP_temp);
    
    iter(j) = mean(numiter(:,1));
    
    BER_SSC(j) = BER_SSC(j)+sum(BER_SSC_temp);
    BLER_SSC(j) = BLER_SSC(j)+sum(BLER_SSC_temp);
end

BER_BP = BER_BP./(N*nbCWs);
BLER_BP = BLER_BP./(nbCWs);

BER_SSC = BER_SSC./(N*nbCWs);
BLER_SSC = BLER_SSC./(nbCWs);
toc;

%% Plots
figure(1);
semilogy(EbN0dB,BER_BP,'-*');hold on;
semilogy(EbN0dB,BLER_BP);
semilogy(EbN0dB,BER_SSC,'-p');
semilogy(EbN0dB,BLER_SSC);
legend('BER - BP', 'FER - BP', 'BER - SSC', 'FER - SSC');
title('BER for bipartite Polar decoding');xlabel('SNRb [dB]');ylabel('BER/FER');
