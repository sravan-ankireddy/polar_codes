#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "functions_SCD_ARM.h"

#include <stdint.h>
#include <inttypes.h>
#include <arm_neon.h>

#define f_macro(L1, L2) sign_macro(L1)*sign_macro(L2)*mini_macro(absl_macro(L1),absl_macro(L2))
#define maxqr 31
#define g_macro(u, L1, L2) mini_macro( maxi_macro( (((1 - 2*u) * L1) + L2), -(maxqr+1)), maxqr)
// #define g_macro(u, L1, L2) (((1 - 2*u) * L1) + L2)
#define sign_macro(x) ((x > 0) - (x < 0))
#define absl_macro(x) (((x > 0) - (x < 0)) * x)
#define mini_macro(x,y) ((x < y) ? x : y)
#define maxi_macro(x,y) ((x < y) ? y : x)

/*is_vec_mem function*/
int is_vec_mem(int in, int Nt, int N, int *data_indices_sorted)
{
    int res = 0;
    int d1 = -1;
    int dN = N/2;
    int i_mem = 0;
    for (i_mem = 0; i_mem < dN; i_mem++)
    {
        if (in == data_indices_sorted[i_mem])
        {
            d1 = i_mem;
            break;
        }
    }

    if (d1>-1)
    {
        int d2;
        d2 = d1 + Nt -1;
        if (d2 < dN && data_indices_sorted[d2] == in+Nt-1)
        {
            res = 1;
        }
    }
    return res; 
}

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

/* swap function int data type */
void swap_int(int* a, int* b) 
{ 
    int t = *a; 
    *a = *b; 
    *b = t; 
}

/* partition function int data type */
int partition_int (int arr[], int low, int high) 
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
        } 
    } 
    swap_int(&arr[i + 1], &arr[high]);
    return (i + 1);
}

/* quicksort function int data type */
void quickSort_int(int arr[], int low, int high) 
{
    if (low < high) 
    { 
        int pi = partition_int(arr, low, high); 
  
        quickSort_int(arr, low, pi - 1); 
        quickSort_int(arr, pi + 1, high); 
    } 
}

/*Encoder*/
void encode(uint8_t *u, int N)
{
    /* No of stages */
    int n = 0;
    /* log function */
    while (N >>= 1) n++;

    uint8x16_t data1, data2;

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
            if (del >= 16)
            {
                for (i_sg = 0; i_sg < del; i_sg+=16)
                {
                    data1 = vld1q_u8(&u[base+i_sg]);
                    data2 = vld1q_u8(&u[base+i_sg+del]);

                    data1 = veorq_u8(data1, data2);

                    vst1q_u8(&u[base+i_sg],data1);
                }
            }

            else
            {
                for (i_sg = 0; i_sg < del; i_sg++)
                {
                    u[base+i_sg] = u[base+i_sg]^u[base+i_sg+del];
                }
            }
        }
    }
}


/* SC Decoder */
void decode(int *msg_cap, int N, int n, int K, int *LLR_Q, int *info_nodes, int *data_pos)
{
    /* Beliefs */
    int L[n+1][N];
    int i_n; 
    int i_N;

    /* Decisions */
    uint8_t ucap[n+1][N];

    /* Node state vector */
    int ns[2*N-1];

    for (i_N = 0; i_N < 2*N-1; i_N++)
    {
        ns[i_N] = 0;
    }

    /* Belief initialisation */
    for (i_N = 0; i_N < N; i_N++)
    {
        L[0][i_N] = LLR_Q[i_N];
    }

    /* Propogation parameters */
    int node = 0;
    int lnode = 0;
    int rnode = 0;
    int depth = 0;
    int ldepth = 0;
    int cdepth = 0;
    int done = 0;
    int npos = 0;
    int temp = 0;
    int ltemp = 0;
    int ctemp = 0;
    int i_L;
    int i_min;
    int node_type_ind;

    int i_s;
    int i_g;
    int del;
    int base;
    int i_sg;
    int n_enc;
    int N_enc;

    /* Traverse till all bits are decoded */
    while (done == 0)
    {
        /* Position of node in node state vector */
        npos = (1 << depth) - 1 + node;

        /* Length of current sub-vector */
        temp = 1 << (n - depth);

        /* Index of current node in node_type vector */
        node_type_ind = node*temp;

        /* Check for leaf node */
        if (depth == n)
        {
            /* Check for frozen node and take decision */
            ucap[n][node] = 0;

            if (info_nodes[node] != 0 && L[n][0] < 0)
                ucap[n][node] = 1;

            (node == N-1) ? (done = 1) :  (node >>= 1 , depth -= 1 );

        }

        /* Non-leaf nodes */
        else
        {
            /* Propogate to left child */
            if (ns[npos] == 0)
            {
                /* f_minsum and storage */
                for (i_L = 0; i_L < temp/2; i_L++)
                {
                    L[depth+1][i_L] = f_macro(L[depth][i_L], L[depth][i_L + temp/2]);
                }

                /* Next node: Left child */
                node <<= 1 ; depth += 1 ;

                /* Incoming belief length for left child */
                temp>>=1;
                
                ns[npos] = 1;
            }
            
            else
            {
                /* Propogate to right child */
                if (ns[npos] == 1)
                {
                    /* g_minsum and storage */
                    for (i_L = 0; i_L < temp/2; i_L++)
                    {
                        L[depth+1][i_L] = g_macro(ucap[depth+1][i_L+node_type_ind], L[depth][i_L ], L[depth][i_L+temp/2]);
                    }

                    /* Next node: right child */
                    node = (node << 1) + 1; depth += 1;

                    /* Incoming belief length for right child */
                    temp >>= 1;

                    ns[npos] = 2;
                }

                /* Propogate to parent node */
                else
                {
                    /* Combine */
                    int count = 0;
                    uint8x16_t data1, data2;
                    if (temp > 16)
                    {
                        for (i_L = 0; i_L < temp/2; i_L+=16)
                        { 
                            data1 = vld1q_u8 (&ucap[depth + 1][i_L + node_type_ind]);
                            data2 = vld1q_u8 (&ucap[depth + 1][i_L + node_type_ind + temp/2]);

                            data1 = veorq_u8 (data1, data2);

                            vst1q_u8(&ucap[depth][i_L + node_type_ind],data1);
                            vst1q_u8(&ucap[depth][i_L + node_type_ind + temp/2],data2);
                        }
                    }

                    else
                    {
                        for (i_L = 0; i_L < temp; i_L++)
                        {
                            if ( count < temp/2)
                            {
                                ucap[depth][i_L + node_type_ind] =  ucap[depth + 1][i_L + node_type_ind]^ucap[depth + 1][i_L + node_type_ind + temp/2];
                            }
                            else
                            {
                                ucap[depth][i_L + node_type_ind] = ucap[depth + 1][i_L + node_type_ind]; 
                            }
                            count++;
                        }
                    }

                    node >>= 1; depth -= 1;

                }
            }
        } // end of non-leaf node else
    } // end of while loop

    int i_m = 0;
    for(i_m = 0; i_m < K; i_m++)
    {
        msg_cap[i_m] = ucap[n][data_pos[i_m]];
    }
}
