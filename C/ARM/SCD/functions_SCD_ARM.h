#include <stdint.h>

int uni ();

float randn (float mu, float sigma);

void swap_int(int* a, int* b);

int partition_int (int arr[], int low, int high);

void quickSort_int(int arr[], int low, int high);

int is_vec_mem(int in, int Nt, int N, int *data_pos_sorted);

void encode(uint8_t* u, int N);

void decode(int *msg_cap, int n, int N, int K, int *LLR_Q, int *info_nodes, int *data_pos);
