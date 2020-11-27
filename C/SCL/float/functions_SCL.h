int uni ();

float randn (float mu, float sigma);

void encode(int *u, int N);

void swap_int(int* a, int* b);

int partition_int (int arr[], int *ind_ord, int low, int high);

void quickSort_int(int arr[], int *ind_ord, int low, int high);

void swap_float(float* a, float* b);

int partition_float (float arr[], int *ind_ord, int low, int high);

void quickSort_float(float arr[], int *ind_ord, int low, int high);

void insertionSort(int arr[], int ind_ord[], int n);

void min_pos(int l, int N, int* arr, int ind_min[]);

void decode_unrolled(int N, int K, int l, int *info_nodes, int *data_pos, float *LLR, int *msg_cap);

void crcGen(int* msg, int K, int* polynomial, int m);

void crcDet(int* msg_cap, int K, int l, int* polynomial, int m, int* crc_check);