%% Simulation Parameters

% Loading Channel indices in increasing order of reliability
Q=[0 1 2 4 8 16 32 3 5 64 9 6 17 10 18 128 12 33 65 20 256 34 24 36 7 129 66 512 11 40 68 130 ...
    19 13 48 14 72 257 21 132 35 258 26 513 80 37 25 22 136 260 264 38 514 96 67 41 144 28 69 42 ...
    516 49 74 272 160 520 288 528 192 544 70 44 131 81 50 73 15 320 133 52 23 134 384 76 137 82 56 27 ...
    97 39 259 84 138 145 261 29 43 98 515 88 140 30 146 71 262 265 161 576 45 100 640 51 148 46 75 266 273 517 104 162 ...
    53 193 152 77 164 768 268 274 518 54 83 57 521 112 135 78 289 194 85 276 522 58 168 139 99 86 60 280 89 290 529 524 ...
    196 141 101 147 176 142 530 321 31 200 90 545 292 322 532 263 149 102 105 304 296 163 92 47 267 385 546 324 208 386 150 153 ...
    165 106 55 328 536 577 548 113 154 79 269 108 578 224 166 519 552 195 270 641 523 275 580 291 59 169 560 114 277 156 87 197 ...
    116 170 61 531 525 642 281 278 526 177 293 388 91 584 769 198 172 120 201 336 62 282 143 103 178 294 93 644 202 592 323 392 ...
    297 770 107 180 151 209 284 648 94 204 298 400 608 352 325 533 155 210 305 547 300 109 184 534 537 115 167 225 326 306 772 157 ...
    656 329 110 117 212 171 776 330 226 549 538 387 308 216 416 271 279 158 337 550 672 118 332 579 540 389 173 121 553 199 784 179 ...
    228 338 312 704 390 174 554 581 393 283 122 448 353 561 203 63 340 394 527 582 556 181 295 285 232 124 205 182 643 562 286 585 ...
    299 354 211 401 185 396 344 586 645 593 535 240 206 95 327 564 800 402 356 307 301 417 213 568 832 588 186 646 404 227 896 594 ...
    418 302 649 771 360 539 111 331 214 309 188 449 217 408 609 596 551 650 229 159 420 310 541 773 610 657 333 119 600 339 218 368 ...
    652 230 391 313 450 542 334 233 555 774 175 123 658 612 341 777 220 314 424 395 673 583 355 287 183 234 125 557 660 616 342 316 ...
    241 778 563 345 452 397 403 207 674 558 785 432 357 187 236 664 624 587 780 705 126 242 565 398 346 456 358 405 303 569 244 595 ...
    189 566 676 361 706 589 215 786 647 348 419 406 464 680 801 362 590 409 570 788 597 572 219 311 708 598 601 651 421 792 802 611 ...
    602 410 231 688 653 248 369 190 364 654 659 335 480 315 221 370 613 422 425 451 614 543 235 412 343 372 775 317 222 426 453 237 ...
    559 833 804 712 834 661 808 779 617 604 433 720 816 836 347 897 243 662 454 318 675 618 898 781 376 428 665 736 567 840 625 238 ...
    359 457 399 787 591 678 434 677 349 245 458 666 620 363 127 191 782 407 436 626 571 465 681 246 707 350 599 668 790 460 249 682 ...
    573 411 803 789 709 365 440 628 689 374 423 466 793 250 371 481 574 413 603 366 468 655 900 805 615 684 710 429 794 252 373 605 ...
    848 690 713 632 482 806 427 904 414 223 663 692 835 619 472 455 796 809 714 721 837 716 864 810 606 912 722 696 377 435 817 319 ...
    621 812 484 430 838 667 488 239 378 459 622 627 437 380 818 461 496 669 679 724 841 629 351 467 438 737 251 462 442 441 469 247 ...
    683 842 738 899 670 783 849 820 728 928 791 367 901 630 685 844 633 711 253 691 824 902 686 740 850 375 444 470 483 415 485 905 ...
    795 473 634 744 852 960 865 693 797 906 715 807 474 636 694 254 717 575 913 798 811 379 697 431 607 489 866 723 486 908 718 813 ...
    476 856 839 725 698 914 752 868 819 814 439 929 490 623 671 739 916 463 843 381 497 930 821 726 961 872 492 631 729 700 443 741 ...
    845 920 382 822 851 730 498 880 742 445 471 635 932 687 903 825 500 846 745 826 732 446 962 936 475 853 867 637 907 487 695 746 ...
    828 753 854 857 504 799 255 964 909 719 477 915 638 748 944 869 491 699 754 858 478 968 383 910 815 976 870 917 727 493 873 701 ...
    931 756 860 499 731 823 922 874 918 502 933 743 760 881 494 702 921 501 876 847 992 447 733 827 934 882 937 963 747 505 855 924 ...
    734 829 965 938 884 506 749 945 966 755 859 940 830 911 871 639 888 479 946 750 969 508 861 757 970 919 875 862 758 948 977 923 ...
    972 761 877 952 495 703 935 978 883 762 503 925 878 735 993 885 939 994 980 926 764 941 967 886 831 947 507 889 984 751 942 996 ...
    971 890 509 949 973 1000 892 950 863 759 1008 510 979 953 763 974 954 879 981 982 927 995 765 956 887 985 997 986 943 891 998 766 ...
    511 988 1001 951 1002 893 975 894 1009 955 1004 1010 957 983 958 987 1012 999 1016 767 989 1003 990 1005 959 1011 1013 895 1006 1014 1017 1018 ...
    991 1020 1007 1015 1019 1021 1022 1023]+1;

% Length of codeword
N = 512;

% Extracting channel indices from Q for snaller N
Q = Q(Q<=N);

% Rate of code
R = 0.5;

% Length of message vector
K_BP_tot = round(N*R);

% Number of blocks transmitted
Nblocks = 10e2;

%% BP parameters
max_iter = 100;
bipartite = 1;

num_par_bits = 8;
num_par_eq = 8;

% K updated
K = K_BP_tot - num_par_eq;
K_BP = K;

A_BP_tot = true(1,N);
A_BP_tot(Q(1:N-K_BP_tot)) = 0;
A_BP_tot = bitrevorder(A_BP_tot);

A_BP = true(1,N);
A_BP(Q(1:N-K_BP)) = 0;
A_BP = bitrevorder(A_BP);

%create H matrix
if bipartite==1
    
    %already kicks out frozen indices
    [H_BP] = polar2bipartite(N,A_BP,1);
    [H_BP_par] = polar2bipartite(N,A_BP_tot,1);
    
    % number of aditional parity bits
    A_info = find(A_BP_tot == 1);
    
    H_par = zeros(num_par_eq,N+num_par_eq);
    
%     del_1 = floor(K_BP_tot/num_par_eq); del_2 = floor(num_par_eq/2);
%     start_ind = 1; end_ind = K_BP_tot;
    
%     par_pos = start_ind:del_1:end_ind;
%     pasr_pos = par_pos(1:num_par_eq);

    del_1 = 20; del_2 = 3;
%     del_1 = 5; del_2 = 1;
    par_pos = 1:del_1:(num_par_eq-1)*del_1 + 1;

    for i = 1:num_par_eq
        H_par(i,par_pos + del_2*(i-1)) = 1;
        H_par(i,N+i) = 1;
    end
    
    %prune all useless nodes
    H_BP_par = pruneGraph(H_BP_par,N);
    H_BP = pruneGraph(H_BP,N);
else
    H_BP_par = createDensePolarH(N,A_BP_tot);
    Hdense = H_BP_par; 
end

[ParRI_BP, ParCI_BP] = find(H_BP == 1);
ParI_BP = [ParRI_BP, ParCI_BP];

decoder_BP = comm.LDPCDecoder('ParityCheckMatrix',ParI_BP,'OutputValue','Whole codeword','DecisionMethod','Soft decision','MaximumIterationCount',max_iter,'NumIterationsOutputPort',1,'IterationTerminationCondition', 'Parity check satisfied');

[ParRI_BP_par, ParCI_BP_par] = find(H_BP_par == 1);
ParI_BP_par = [ParRI_BP_par, ParCI_BP_par];

decoder_BP_par = comm.LDPCDecoder('ParityCheckMatrix',ParI_BP_par,'OutputValue','Whole codeword','DecisionMethod','Soft decision','MaximumIterationCount',max_iter,'NumIterationsOutputPort',1,'IterationTerminationCondition', 'Parity check satisfied');

H_tot = zeros(size(H_BP_par,1)+num_par_eq, size(H_BP_par,2)+num_par_eq);
H_tot(1:size(H_BP_par,1),1:size(H_BP_par,2)) = H_BP_par;
H_tot(size(H_BP_par,1)+1:end,end-N-num_par_eq+1:end) = H_par;

[ParRI_tot, ParCI_tot] = find(H_tot == 1);
ParI_tot = [ParRI_tot, ParCI_tot];

decoder_tot = comm.LDPCDecoder('ParityCheckMatrix',ParI_tot,'OutputValue','Whole codeword','DecisionMethod','Soft decision','MaximumIterationCount',max_iter,'NumIterationsOutputPort',1,'IterationTerminationCondition', 'Parity check satisfied');

%%

% Energy per bit
Eb = 1;

% No. of samples in EBN0 range
nSam = 8;

% start and stop of EbN0dB
EbN0dB_start = 1; EbN0dB_stop = 3.7;

% Eb/N0 step size
delta = (EbN0dB_stop - EbN0dB_start)/(nSam-1);

% Eb/N0 range of signal in dB
EbN0dB = EbN0dB_start:delta:EbN0dB_stop;

% Varinace of noise in linear scale
var_N = (10.^(-EbN0dB/10))/(2 * R);

% Vector to store Bit Error Rate
BER_SSCD = zeros(1,nSam);
BER_BP = zeros(1,nSam);
BER_BP_tot = zeros(1,nSam);

% Vector to store Bit Error Rate
BLER_SSCD = zeros(1,nSam);
BLER_BP = zeros(1,nSam);
BLER_BP_tot = zeros(1,nSam);

% Indices corresponding to data channels
data_pos = Q(N-K+1:end);

% Boolean array with information nodes pos = 1
info_check_vec = zeros(1,N);
info_check_vec(data_pos) = 1;

test = zeros(log2(N)+1,N);
node_type_mat = find_node_type(test, info_check_vec, 1, N);

tic;
% Simulation over the range of EbN0
for i_E = 1:length(EbN0dB)

	% Display the simulation status
	disp("counter " + i_E + " of " + length(EbN0dB));
    disp(' ');

	% Std Dev of AWGN noise
	sig = sqrt(var_N(i_E));
    
    %Initialising the number of errors to zero for each sigma
    err_SSCD = 0;
    fer_SSCD = 0;
    err_BP = 0;
    fer_BP = 0;
    err_BP_tot = 0;
    fer_BP_tot = 0;
    
    Nblocks_count = 0;
	% Simulation over the message vectors
	parfor i_m = 1:Nblocks
        
        Nblocks_count = Nblocks_count + 1;
		% Generating a random message vector
        msg_BP_tot = randi([0,1],K_BP_tot,1);
        msg_SSC = msg_BP_tot(1:K);
        msg_BP = msg_BP_tot(1:K_BP);

		% Tx vector initialisation
		u_SSC = zeros(N,1);
        u_BP_tot = zeros(N,1);
        u_BP = zeros(N,1);

		% Assigning data to data channels
		u_SSC(data_pos) = msg_SSC;
        u_BP_tot(A_BP_tot) = msg_BP_tot;
        u_BP(A_BP) = msg_BP;
               
		% Encoding the Tx vector
		x_SSC = encode(u_SSC);
        x_BP = encode(u_BP);
        x_BP_tot = encode(u_BP_tot);
        
        % add parity
        x_BP_tot = [x_BP_tot; zeros(num_par_eq,1)];
        for i_par = 1:num_par_eq
            x_BP_tot(N+i_par) = mod(sum(x_BP_tot(par_pos+ del_2*(i_par-1))),2);
        end
        
        % BPSK modulation
		x_SSC = sqrt(Eb)*(-1).^x_SSC;
        x_BP = (1 - 2*x_BP);
        x_BP_tot = (1 - 2*x_BP_tot);

		% Adding AWGN noise
        noise_BP_par = sig * randn(size(x_BP_tot));
        y_BP_tot = x_BP_tot + noise_BP_par;
        
        noise_SSC = noise_BP_par(1:N);
		y_SSC = x_SSC + noise_SSC;
        
        noise_BP = noise_BP_par(1:N);
        y_BP = x_BP + noise_BP;
        
        LLR_SSC = 2*y_SSC/(sig^2);
        LLR_BP_tot = 2*y_BP_tot/(sig^2);
        LLR_BP = 2*y_BP/(sig^2);
        
        % SC Decoding
        msg_hat_SSC = decode_SSCD(LLR_SSC, N, node_type_mat, info_check_vec, data_pos);
        
        % BP decoding
        Lch_ext_BP = zeros(size(H_BP,2),1);
        Lch_ext_tot = zeros(size(H_tot,2),1);
        
        % Last N correspond to received vector
        Lch_ext_BP((end-N+1):end) = LLR_BP;
        LLR_BP_hat = step(decoder_BP, Lch_ext_BP); xhat_BP=(-sign(LLR_BP_hat))/2+0.5;
        
        Lch_ext_tot((end-N-num_par_eq+1):end) = LLR_BP_tot;
        LLR_tot_hat = step(decoder_tot, Lch_ext_tot); xhat_tot=(-sign(LLR_tot_hat))/2+0.5;
                
        msg_hat_BP = xhat_BP(end-N+1:end); msg_hat_BP = polarTransform(msg_hat_BP, A_BP); msg_hat_BP = msg_hat_BP(A_BP);
        msg_hat_tot = xhat_tot(end-N-num_par_eq+1:end-num_par_eq); msg_hat_tot = polarTransform(msg_hat_tot, A_BP_tot); msg_hat_tot = msg_hat_tot(A_BP_tot);
        
        % Checking for errors
        nberr_SSC = sum(xor(msg_SSC,msg_hat_SSC));
        err_SSCD = err_SSCD + nberr_SSC;
        if (nberr_SSC > 0)
            fer_SSCD = fer_SSCD + 1;
        end
        
        nberr_BP = sum(xor(msg_BP,msg_hat_BP));
        err_BP= err_BP + nberr_BP;
        if (nberr_BP > 0)
            fer_BP = fer_BP + 1;
        end
        
        nberr_BP_tot = sum(xor(msg_BP_tot,msg_hat_tot));
        err_BP_tot = err_BP_tot + nberr_BP_tot;
        if (nberr_BP_tot > 0)
            fer_BP_tot = fer_BP_tot + 1;
        end
                
    end    
    % Bit Error Rate & Frame Error Rate
    BER_SSCD(1,i_E) = err_SSCD/(Nblocks_count*N);
    BLER_SSCD(1,i_E) = fer_SSCD/(Nblocks_count);
    
    BER_BP(1,i_E) = err_BP/(Nblocks_count*N);
    BLER_BP(1,i_E) = fer_BP/(Nblocks_count);
    
    BER_BP_tot(1,i_E) = err_BP_tot/(Nblocks_count*N);
    BLER_BP_tot(1,i_E) = fer_BP_tot/(Nblocks_count);

end
time = toc;
disp("Time taken for length " + N + " and " + Nblocks*nSam + " simulations is " + time + " secs");
%% plot for BER_SSCL vs Eb/N0
figure(1)
semilogy(EbN0dB,BER_SSCD,'-s');
hold on;
semilogy(EbN0dB,BLER_SSCD,'-s');
semilogy(EbN0dB,BER_BP,'--*');
semilogy(EbN0dB,BLER_BP,'--*');
semilogy(EbN0dB,BER_BP_tot,'--x');
semilogy(EbN0dB,BLER_BP_tot,'--x');

legend('BER SSC','BLER SSC','BER BP','BLER BP', 'BER BP par','BLER BP par');
str = sprintf('Plot for BER vs Eb/N0 for SSCD and BP');
title(str);
xlabel('Eb/N0 in dB');
ylabel('Bit Error Rate(BER)');