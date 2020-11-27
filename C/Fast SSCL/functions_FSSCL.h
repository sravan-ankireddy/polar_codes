int is_vec_mem(int in, int Nt, int N, int *data_pos_sorted);

int uni ();

float randn (float mu, float sigma);

void encode(int *u, int N);

void swap_int(int* a, int* b);

int partition_int (int arr[], int *ind_ord, int low, int high);

void quickSort_int(int arr[], int *ind_ord, int low, int high);

void insertionSort(int arr[], int ind_ord[], int n);

void min_pos(int l, int N, int* arr, int ind_min[]);

void find_node_type(int *node_type, int N, int Nt, int in, int *data_pos_sorted, int *frozen_pos_sorted);

void decode_unrolled(int N, int K, int l, int *info_nodes, int *data_pos, int *node_type, int *LLR_Q, int *msg_cap);

void crcGen(int* msg, int K, int* polynomial, int m);

void crcDet(int* msg_cap, int K, int l, int* polynomial, int m, int* crc_check);

void dec2bin(int i, int *str, int n);