#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "functions_SCL.h"
#define maxqr 31
#define sign_macro(x) ((x > 0) - (x < 0))
#define f_macro(L1, L2) sign_macro(L1)*sign_macro(L2)*mini_macro(absl_macro(L1),absl_macro(L2))
// #define g_macro(u, L1, L2) mini_macro( maxi_macro( (((1 - 2*u) * L1) + L2), -(maxqr+1)), maxqr)
#define g_macro(u, L1, L2) (((1 - 2*u) * L1) + L2)
#define absl_macro(x) (((x > 0) - (x < 0)) * x)
#define mini_macro(x,y) ((x < y) ? x : y)
#define maxi_macro(x,y) ((x < y) ? y : x)

/*Uniform distribution number generator*/
int uni ()
{
    return (rand()%2);
}

/*Gaussian random number generator*/
float randn (float mu, float sigma)
{
    float U1, U2, W, mult;
    static float X1, X2;
    static int call = 0;

    if (call == 1)
    {
      call = !call;
      return (mu + sigma * (float) X2);
    }

    do
    {
      U1 = -1 + ((float) rand () / RAND_MAX) * 2;
      U2 = -1 + ((float) rand () / RAND_MAX) * 2;
      W = pow (U1, 2) + pow (U2, 2);
    }
    while (W >= 1 || W == 0);

    mult = sqrt ((-2 * log (W)) / W);
    X1 = U1 * mult;
    X2 = U2 * mult;

    call = !call;

    return (mu + sigma * (float) X1);
}

/*Encoder*/
void encode(int *u, int N)
{
    /* No of stages */
    int n = 0;
    /* log function */
    while (N >>= 1) n++;

    int i_s;
    int i_g;
    int del;
    int base;
    int i_sg;

    /* Stage after stage */
    for (i_s = 0; i_s < n; i_s++)
    {
        del = (1 << i_s);

        /* Group after group in each stage */
        for (i_g = 0; i_g < (1 << (n-i_s-1)); i_g++)
        {
            base = (1 << (i_s+1)) * (i_g);

            /* Sub group after sub group */
            for (i_sg = 0; i_sg < del; i_sg++)
            {
                u[base+i_sg] = u[base+i_sg]^u[base+i_sg+del];
            }
        }
    }
}

/* SCL Decoder */
void decode_unrolled(int N, int K, int l, int *info_nodes, int *data_pos, float *LLR_Q, int *msg_cap)
{
    int n = 0, N_log = N, l0_log = l, ln = 0;
     /* log function */
    while (N_log >>= 1) n++;
    while (l0_log >>= 1) ln++;

    /* Beliefs */
    float L[n+1][l*N];
    /* Decisions */
    int beta[n+1][l*N];
    /* Orders of surviving decoders */
    int ind_ord_mat[n+1][l*N];

    /* Propagation Parameters */
    int i_list;
    int i_L;
    int i_N;
    int node = 0;
    int depth = 0;
    int done = 0;
    int npos = 0;
    int temp = 0;
    int node_type_ind;
    
    int ind_ord[2*l];
    int ind_ord_temp[2*l];
    int ind_ord_l[l];
    int codeword[l];
    int codeword_temp[2*l];
    float PM[l];
    float PM_temp[2*l];

    int start_ind, end_ind;
    int i_temp;
    int del;
    int base;
    int i_s;
    int i_g;
    int i_sg;

    float key_arr;
    int key_ind, cur_index;

    /* Node state vector */
    int ns[2*N-1];

    for (i_N = 0; i_N < 2*N-1; i_N++)
    {
        ns[i_N] = 0;
    }

    /* Belief initialization */
    for (i_list = 0; i_list < l; i_list++)
    {
        for (i_N = 0; i_N < N; i_N++)
        {
            L[0][i_N + i_list*N] = LLR_Q[i_N];
        }
    }

    /* Initializing Path Metrics and order variable */
    for (i_list = 0; i_list < l; i_list++)
    {
        PM[i_list] = 0;
    }

    /* Variable to count no. of decoded data bits */
    int counter = 0;

    while (!(done == 1 && depth == -1)) //traverse till all bits are decoded and root node is reached again
    {
        /* Position of node in node state vector */
        npos = (1 << depth) - 1 + node;

        /* Length of current sub-vector */
        temp = 1 << (n - depth);

        /* Index of current node in node_type vector */
        node_type_ind = node*temp;

        /* Check for leaf node */
        if (depth == n && done == 0)
        {
            if ( info_nodes[node_type_ind] == 0 )
            {
                /* Assigning 0 to decoded bits and updating path metrics */
                for (i_list = 0; i_list < l; i_list++)
                {
                    beta[depth][node_type_ind + i_list*N] = 0;
                    
                    ( L[depth][i_list*N] < 0 ) ? ( PM[i_list] = PM[i_list] - L[depth][i_list*N] ) : 1 ;
                }
            
                for (i_list = 0; i_list < l; i_list++)
                {
                    ind_ord_l[i_list] = i_list;
                }

                /* Sorting the Path Metrics -- Insertion Sort */
                for (i_list = 1; i_list < l; i_list++) 
                {
                    key_arr = PM[i_list];
                    key_ind = ind_ord_l[i_list]; 
                    cur_index = i_list-1; 
 
                    while (cur_index >= 0 && PM[cur_index] > key_arr) 
                    { 
                        PM[cur_index+1] = PM[cur_index];
                        ind_ord_l[cur_index+1] = ind_ord_l[cur_index];
                        cur_index = cur_index-1; 
                    } 
                    PM[cur_index+1] = key_arr;
                    ind_ord_l[cur_index+1] = key_ind;
                }

                /* New ordering */
                for (i_list = 0; i_list < l; i_list++)
                {
                    ind_ord_mat[depth][node + i_list*N] = ind_ord_l[i_list];
                }

            } // end of frozen leaf node check

            else
            {
                if (counter > ln - 1)
                {
                    /* Duplicating the old data */
                    for (i_list = 0; i_list < l; i_list++)
                    {
                        PM_temp[i_list] = PM[i_list];
                        PM_temp[i_list + l] = PM[i_list];
                    }

                    /* Forking into Path 0 and Path 1 */
                    for (i_list = 0; i_list < l; i_list++)
                    {
                        codeword_temp[i_list] = 0;
                        
                        if ( L[depth][i_list*N] < 0 )
                        {
                            PM_temp[i_list] = PM_temp[i_list] - L[depth][i_list*N];
                        }
                    }

                    for (i_list = l; i_list < 2*l; i_list++)
                    {
                        codeword_temp[i_list] = 1;
                        
                        if (  L[depth][(i_list-l)*N] > 0 )
                        {
                            PM_temp[i_list] = PM_temp[i_list] + L[depth][(i_list-l)*N];
                        }
                    }

                    
                    for (i_list = 0; i_list < 2*l; i_list++)
                    {
                        ind_ord[i_list] = i_list;
                    }

                    /* Sorting the Path Metrics -- Insertion Sort */
                    for (i_list = 1; i_list < 2*l; i_list++) 
                    {

                        key_arr = PM_temp[i_list];
                        key_ind = ind_ord[i_list]; 
                        cur_index = i_list-1; 
     
                        while (cur_index >= 0 && PM_temp[cur_index] > key_arr) 
                        { 
                            PM_temp[cur_index+1] = PM_temp[cur_index];
                            ind_ord[cur_index+1] = ind_ord[cur_index];
                            cur_index = cur_index-1; 
                        } 
                        PM_temp[cur_index+1] = key_arr;
                        ind_ord[cur_index+1] = key_ind;
                    }

                    for (i_list = 0; i_list < l; i_list++)
                    {
                        PM[i_list] = PM_temp[i_list];
                    }
                    
                    /* Re-ordering the estimated vector and updating beta*/
                    for (i_list = 0; i_list < l; i_list++)
                    {
                        beta[depth][node_type_ind + i_list*N] = codeword_temp[ind_ord[i_list]];

                        ind_ord_mat[depth][node + i_list*N] = ind_ord[i_list]%l;
                    }

                    counter++;
                } // end of normal data leaf node check

                else
                {
                    int init;

                    if (counter > 0)
                    {
                        for (init = 0; init < 2*counter; init++)
                        {
                            start_ind = 0 + init*(l/(2*counter)); end_ind = start_ind + l/(4*counter);

                            for (i_list = start_ind; i_list < end_ind; i_list++)
                            {
                                codeword[i_list] = 0;
                            }
                            start_ind = end_ind; end_ind = start_ind + l/(4*counter);

                            for (i_list = start_ind; i_list < end_ind; i_list++)
                            {
                                codeword[i_list] = 1;
                            }
                        }
                    }

                    else
                    {
                        for (i_list = 0; i_list < l/2; i_list++)
                        {
                            codeword[i_list] = 0;
                        }
                        for (i_list = l/2; i_list < l; i_list++)
                        {
                            codeword[i_list] = 1;
                        }
                    }

                    /* Updating the Path Metric */
                    for (i_list = 0; i_list < l; i_list++)
                    {
                        if (codeword[i_list] != 0.5*(1 - sign_macro(L[depth][i_list*N])))
                        {
                            PM[i_list] = PM[i_list] + absl_macro(L[depth][i_list*N]);
                        }
                    }

                    /* Initialization for new estimate */
                    for (i_list = 0; i_list < l; i_list++)
                    {
                        ind_ord_l[i_list] = i_list;
                    }

                    /* Sorting the Path Metrics  */
                    for (i_list = 1; i_list < l; i_list++) 
                    {

                        key_arr = PM[i_list];
                        key_ind = ind_ord_l[i_list]; 
                        cur_index = i_list-1; 
     
                        while (cur_index >= 0 && PM[cur_index] > key_arr) 
                        { 
                            PM[cur_index+1] = PM[cur_index];
                            ind_ord_l[cur_index+1] = ind_ord_l[cur_index];
                            cur_index = cur_index-1; 
                        } 
                        PM[cur_index+1] = key_arr;
                        ind_ord_l[cur_index+1] = key_ind;
                    }

                    /* New ordering */
                    for (i_list = 0; i_list < l; i_list++)
                    {
                        ind_ord[i_list] = ind_ord_l[i_list];
                        ind_ord[i_list + l] = ind_ord_l[i_list];

                        ind_ord_mat[depth][node + i_list*N] = ind_ord_l[i_list];
                    }

                    /* Re-ordering the estimated vector */
                    for (i_list = 0; i_list < l; i_list++)
                    {
                        beta[depth][node_type_ind + i_list*N] = codeword[ind_ord[i_list]];
                    }

                    counter++;
                } //end of data leaf node check

            } // end of else 

            (node == N-1) ? (done = 1, node >>= 1 , depth -= 1) :  (node >>= 1 , depth -= 1 );
        } // end of leaf node

        else
        {
            /* Propogate to left child */
            if (ns[npos] == 0 && done == 0)
            {
                int i_temp;

                /* f_minsum and storage */
                for (i_list = 0; i_list < l; i_list++)
                {
                    i_temp = i_list*N;

                    for (i_L = 0; i_L < temp/2; i_L++)
                    {
                        L[depth + 1][i_L + i_temp] = f_macro(L[depth][i_L + i_temp], L[depth][i_L + temp/2 + i_temp]);
                    }
                }

                /* Next node: Left child */
                node <<= 1 ; depth += 1 ;

                /* Incoming belief length for left child */
                temp>>=1;

                ns[npos] = 1;
            } // end of left child propogation

            else
            {
                /* Propogate to right child */               
                if (ns[npos] == 1 && done == 0)
                {
                    int i_temp;
                    int i_temp_ord;

                    /* g_minsum and storage */
                    for (i_list = 0; i_list < l; i_list++)
                    {
                        i_temp = i_list*N;
                        i_temp_ord = (ind_ord_mat[depth+1][2*node + i_list*N])*N;

                        for (i_L = 0; i_L < temp/2; i_L++)
                        {
                            L[depth + 1][i_L + i_temp] = g_macro(beta[depth + 1][i_L + node_type_ind + i_temp], L[depth][i_L + i_temp_ord], L[depth][i_L + temp/2 + i_temp_ord]);
                        }
                    }

                    /* Next node: right child */
                    node = (node << 1) + 1; depth += 1;

                    /* Incoming belief length for right child */
                    temp >>= 1;

                    ns[npos] = 2;
                } // end of right child propogation

                /* Propogate to parent node */
                else
                {
                    /* Updating beta */
                    int i_up;
                    int i_temp_ord;
                    int i_temp;
                    int ind_ord_temp2[l];

                    for (i_list = 0; i_list < l; i_list++)
                    {
                        /* Order of surviving decoders from left */
                        ind_ord_temp[i_list] = ind_ord_mat[depth +1][2*node + i_list*N];

                        /* Order of surviving decoders from right */
                        ind_ord_temp2[i_list] = ind_ord_mat[depth + 1][2*node + 1 + i_list*N];
                    }

                    /* Final Order of surviving decoders to be passed to parent node */
                    for (i_list = 0; i_list < l; i_list++)
                    {
                        i_temp = ind_ord_temp2[i_list];

                        ind_ord_mat[depth][node + i_list*N] = ind_ord_temp[i_temp];
                    }

                    /* Updating beta to be passed to parent node */
                    for (i_list = 0; i_list < l; i_list++)
                    {
                        i_temp = node_type_ind + i_list*N;
                        i_temp_ord = node_type_ind + (ind_ord_mat[depth + 1][2*node + 1 + i_list*N])*N;

                        for (i_up = 0; i_up < temp/2; i_up++)
                        {
                            beta[depth][i_up + i_temp] = beta[depth+1][i_up + i_temp_ord] ^ beta[depth+1][i_temp + temp/2 + i_up];
                        }

                        for (i_up = 0; i_up < temp/2; i_up++)
                        {
                            beta[depth][i_up + temp/2 + i_temp] =  beta[depth+1][i_temp + temp/2 + i_up];
                        }
                    }

                    node >>= 1; depth -= 1;
                } // end of parent node propagation

            } // end of else for left child propagation

        } // end of non-leaf node check

    }   // end of while loop ==> decoder stops

    int temp_beta[N];
    int i_n;   
    int i_m = 0;

    for (i_list = 0; i_list < l; i_list++)
    {
        i_n = i_list*N;
        for (i_N = 0; i_N < N; i_N++)
        {
            temp_beta[i_N] = beta[0][i_n + i_N];
        }

        /* encoding the decidions to get information at corresponding leaf nodes*/
        for (i_s = 0; i_s < n; i_s++)
        {
            del = (1 << i_s);

            /* Group after group in each stage */
            for (i_g = 0; i_g < (1 << (n-i_s-1)); i_g++)
            {
                base = (del << 1)*i_g;

                /* Sub group after sub group */
                for (i_sg = 0; i_sg < del; i_sg++)
                {
                    temp_beta[base + i_sg] = temp_beta[base + i_sg]^temp_beta[base + i_sg + del];
                }
            }
        }

        i_n = i_list*K;
        for(i_m = 0; i_m < K; i_m++)
        {
            msg_cap[i_n + i_m] = temp_beta[data_pos[i_m]];
        }
    }
}

/* quickSort_int function for int data type */
void quickSort_int(int arr[], int ind_ord[], int low, int high) 
{
    if (low < high) 
    { 
        int pi = partition_int(arr, ind_ord, low, high); 

        quickSort_int(arr, ind_ord, low, pi - 1); 
        quickSort_int(arr, ind_ord, pi + 1, high); 
    } 
}

/* quickSort_int function for int data type */
void quickSort_float(float arr[], int ind_ord[], int low, int high) 
{
    if (low < high) 
    { 
        int pi = partition_float(arr, ind_ord, low, high); 

        quickSort_float(arr, ind_ord, low, pi - 1); 
        quickSort_float(arr, ind_ord, pi + 1, high); 
    } 
}

/* swap function for integer data type */
void swap_int(int* a, int* b) 
{ 
    int t = *a; 
    *a = *b; 
    *b = t; 
}

/* swap function for integer data type */
void swap_float(float* a, float* b) 
{ 
    float t = *a; 
    *a = *b; 
    *b = t; 
}

/* Partition function for integer data type */
int partition_int(int arr[], int ind_ord[], int low, int high) 
{ 
    int pivot = arr[high]; 
    int i = (low - 1);
  
    int j = low;
    for (j = low; j <= high - 1; j++) 
    { 
        if (arr[j] <= pivot) 
        { 
            i++; 
            swap_int(&arr[i], &arr[j]);
            swap_int(&ind_ord[i], &ind_ord[j]);
        } 
    } 
    swap_int(&arr[i + 1], &arr[high]);
    swap_int(&ind_ord[i + 1], &ind_ord[high]);
    return (i + 1);
}

/* Partition function for integer data type */
int partition_float(float arr[], int ind_ord[], int low, int high) 
{ 
    float pivot = arr[high]; 
    int i = (low - 1);
  
    int j = low;
    for (j = low; j <= high - 1; j++) 
    { 
        if (arr[j] <= pivot) 
        { 
            i++; 
            swap_float(&arr[i], &arr[j]);
            swap_int(&ind_ord[i], &ind_ord[j]);
        } 
    } 
    swap_float(&arr[i + 1], &arr[high]);
    swap_int(&ind_ord[i + 1], &ind_ord[high]);
    return (i + 1);
}

/* Function to find position of minimum element in a vector */
void min_pos(int l, int N, int* arr, int ind_min[])
{
    int min_val;
    int i_list, i_N;
    int temp = N/l;

    for (i_list = 0; i_list < l; i_list++)
    {
        ind_min[i_list] = 0;
        min_val = absl_macro(arr[0 + i_list*temp]);
        for (i_N = 1; i_N < temp; i_N++)
        {
            if (absl_macro(arr[i_N + i_list*temp]) < min_val)
            {
                min_val = absl_macro(arr[i_N + i_list*temp]);
                ind_min[i_list] = i_N;
            }
        }
    }

}

void crcGen(int* msg, int K, int* polynomial, int m)
{
    int i;
    int msg_copy[K];
    for (i = 0; i < K; i++)
    {
        msg_copy[i] = msg[i];
    }

    int j;
    for ( i = 0; i < K - m + 1; i++)
    {
        if ( msg_copy[i] != 0)
        {
            for ( j = 0; j <  m; j++)
            {
                msg_copy[i + j] = msg_copy[i + j] ^ polynomial[j];
            }
        }
    }

    for ( i = K - m; i < K; i++)
    {
        msg[i] = msg_copy[i];
    }
}

void crcDet(int* msg_cap, int K, int l, int* polynomial, int m, int* crc_check)
{
    int il,i,j;
    int sum;
    int msg_copy[K];

    for(il = 0; il < l; il++)
    {
        crc_check[il] = 0;

        for (i = 0; i < K; i++)
        {
            msg_copy[i] = msg_cap[i + il*K];
        }
        
        for ( i = 0; i < K - m + 1; i++)
        {
            if ( msg_copy[i] != 0)
            {
                for ( j = 0; j <  m; j++)
                {
                    msg_copy[i + j] = msg_copy[i + j] ^ polynomial[j];
                }
            }
        }

        sum = 0;
        for ( i = K - m; i < K; i++)
        {
            sum += msg_copy[i];
        }

        if (sum > 0)
        {
            crc_check[il] = 1;
        }
    }
}

/* Function to sort an array using insertion sort*/
void insertionSort(int arr[], int ind_ord[], int n) 
{ 
   int i_in, key_arr, key_ind, cur_index; 
   for (i_in = 1; i_in < n; i_in++) 
   { 
       key_arr = arr[i_in];
       key_ind = ind_ord[i_in]; 
       cur_index = i_in-1; 

       while (cur_index >= 0 && arr[cur_index] > key_arr) 
       { 
           arr[cur_index+1] = arr[cur_index];
           ind_ord[cur_index+1] = ind_ord[cur_index];
           cur_index = cur_index-1; 
       } 
       arr[cur_index+1] = key_arr;
       ind_ord[cur_index+1] = key_ind;
   } 
}
