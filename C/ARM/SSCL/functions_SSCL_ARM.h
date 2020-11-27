#include <stdint.h>

int is_vec_mem(int in, int Nt, int N, int *data_pos_sorted);

int uni ();

float randn (float mu, float sigma);

void encode(uint8_t *u, int N);

void swap_int(int* a, int* b);

int partition_int (int arr[], int *ind_ord, int low, int high);

void quickSort_int(int arr[], int *ind_ord, int low, int high);

void insertionSort(int arr[], int ind_ord[], int n);

void min_pos(int l, int N, int* arr, int ind_min[]);

void find_node_type(int *node_type, int N, int Nt, int in, int *data_pos_sorted, int *frozen_pos_sorted);

void decode_unrolled(int N, int K, int l, int *info_nodes, int *data_pos, int *node_type, int *LLR_Q, uint8_t *msg_cap);

void crcGen(uint8_t* msg, int K, uint8_t* polynomial, int m);

void crcDet(uint8_t* msg_cap, int K, int l, uint8_t* polynomial, int m, int* crc_check);