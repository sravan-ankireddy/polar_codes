# polar_codes
Construction, encoding and decoding of Polar Codes. List decoding of Polar Codes.

This project aims to develop efficient encoder and decoder algorithms for polar codes. In particular, we studied the two limitations of decoders for polar codes: error correction performance and computational complexity. The first decoding algorithm proposed for polar codes, successive cancellation (SC) decoding, performs very poorly compared to its LDPC counterparts. To improve this, we implemented a list decoding algorithm, a generalization of SC decoder where L most probable decoding paths are considered concurrently at each decoding stage. Additionally, a short length Cyclic Redundancy Check (CRC) is added to identify the correct codeword among L estimates. To reduce the computational complexity, a number of simplifications based on the identification of special nodes have been incorporated, and the Log-Likelihood Ratio (LLR) update functions have been optimized to reduce the amount of data copying, resulting in improved decoding speed. Finally, we explored some novel ideas like concatenating Low-Density Parity Check (LDPC) codes with Polar Codes to achieve superior decoding performance with a trade-off for computational complexity. And using a modified Belief Propagation (BP) decoder to traverse the Tanner graph structure of Polar Codes, which enables the reuse of optimized architecture of Message Passing Algorithm (MPA) with a trade-off for error-correcting performance.

Key words: Polar Codes; Successive Cancellation Decoding; List Decoding; Belief Propagation; 5G NR

Key instructions/commands for C code:

Instructions to run C codes on x86 processor:

Change the extensions for functions_ xxxx.c  and simulate_xxxx.c according to the decoder. This an example for SCL.

=======> gcc -O3 functions_SCL.c -o hello simulate_SCL.c -lm

Instructions to run C codes on ARM processor:

Change the extensions for functions_ and simulate_ according to the decoder. This an example for SCL.

=======> gcc -O3 -march=armv8-a -mtune=cortex-a53 -mfpu=neon -ftree-vectorize functions_SCL.c -o hello simulate_SCL.c -lm

=====================================================

functions_SCL.h — header file

functions_SCL.c — all the functions are defined in this file

simulate_SCL.c — all simulation parameters are defined in this function

Abbrevations : 

AWGN - Additive White Gaussian Noise

BCH - Bose-Chaudhuri-Hocquenghem

BEC - Binary Erasure Channel

BER - Bit Error Rate

BLER - Block Error Rate

BPSK - Binary Phase Shift Keying

BSC - Binary Symmetric Channel

CRC - Cyclic Redundancy Check

CRC-BP - CRC aided Belief Propagation Decoder

FER - Frame Error Rate

FFT - Fast Fourier Transform

GF - Galio Field

HWF - Hardware Friendly

LDPC - Low Density Parity Check

LLR - Log Likelihood Ratio

PM - Path Metric

Rep - Repetition

RS - Reed-Solomon

SC - Successive Cancellation

SCL - Successive Cancellation List

SNR - Signal to Noise Ratio

SPC - Single Parity Check

SSC - Simplified Successive Cancellation

SSCL - Simplified Successive Cancellation List
