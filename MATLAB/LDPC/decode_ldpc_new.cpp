
/*LDPC Decoder*/

#include "mex.h"
#include "matrix.h"			// for Matlab mx and mex fuctions
#include "math.h"
#include <stdlib.h>			// what for 
#include "decodeutil_new.h"

#define INF 1000  //maximum value of check_node to variable_node LLR 


double atanh(double x);


//decode(max_iter,vhat,mrows,ncols,iter_out,gamma_n,check_node_ones,max_check_degree,BIGVALUE_COLS,
//                                variable_node_ones,max_variable_degree,BIGVALUE_ROWS);

void decode(double max_iter, double *vhat, int mrows, int ncols, double *iter_out, double *gamma_n,
			double *check_node_ones,double max_check_degree, double BIGVALUE_COLS,
			double *variable_node_ones,double max_variable_degree,double BIGVALUE_ROWS)

{//function braces
	
	int i=0, j=0;
	
	double **sg_array,**sa_array, **bitmessage_temp ;
	double **H, **check_node_ones_matrix,**variable_node_ones_matrix;
	
	sg_array		= matrix(0,mrows-1,0,ncols-1);	
	sa_array		= matrix(0,mrows-1,0,ncols-1);	
	H				= matrix(0,mrows-1,0,ncols-1);	

	check_node_ones_matrix=matrix(0,mrows-1,0,max_check_degree-1);	//for matrix from Matlab containing for each row the column indeces which are '1'
	variable_node_ones_matrix=matrix(0,max_variable_degree-1,0,ncols-1);



		for (i=0;i<max_check_degree;i++)
		{
			for(j=0;j<mrows;j++)
			{
				check_node_ones_matrix[j][i]=*(check_node_ones++); //writing out check_node_ones in 2D matrix form, from Matlab it is passed as one long vector
			
			}
		}

		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		// to print to screen the 'check_node_ones_matrix' matrix
		/*	for ( i=0;i<mrows;i++)
		{
			mexPrintf("\n");
			for( j=0;j<max_check_degree;j++)
			{
				mexPrintf("\t%f",check_node_ones_matrix[i][j]);
			}
			mexPrintf("\n");
		} */
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


		for(i=0;i<ncols;i++)
		{
			for(j=0;j<max_variable_degree;j++)
			{
				variable_node_ones_matrix[j][i]=*(variable_node_ones++);
			}
			
		}


		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		// to print to screen the 'variable_node_ones_matrix' matrix
	/*	for ( i=0;i<max_variable_degree;i++)
		{
			mexPrintf("\n");
			for( j=0;j<ncols;j++)
			{
				mexPrintf("\t%f",variable_node_ones_matrix[i][j]);
			}
			mexPrintf("\n");
		} */
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		
		
		double element=0,col_index=0,row_index=0;
		int temp_row_index=0;

		//initializing the matrices
		for ( i=0;i<ncols;i++)
		{
			for(int j=0;j<max_variable_degree;j++)
			{
				
				row_index=variable_node_ones_matrix[j][i];
				if(row_index==BIGVALUE_ROWS)
					break;


				temp_row_index=(int)row_index;

				sg_array[temp_row_index][i]=gamma_n[i];
				sa_array[temp_row_index][i]=0;
				H[temp_row_index][i]=1;
	
					//mexPrintf("\t%f",sg_array[j][i]);
			}
		}



		//old way of initializing the matrices
/*		for ( i=0;i<ncols;i++)
		{
			for(j=0;j<mrows;j++)
			{
				sg_array[j][i]=*(sg++);   //writing out sg_array in 2D matrix form, from Matlab it is passed as one long vector
				sa_array[j][i]=0;
				H[j][i]=*(h++);				//writing out H in 2D matrix form, from Matlab it is passed as one long vector
			}
		}*/





		

			
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		// to print to screen the 'H' matrix
	/*	for ( i=0;i<mrows;i++)
		{
			mexPrintf("\n");
			for(int j=0;j<ncols;j++)
			{
				mexPrintf("\t%f",H[i][j]);
			}
			mexPrintf("\n");
		} */
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




	bitmessage_temp = matrix(0,mrows-1,0,ncols-1);
	double *sum_of_b;
	sum_of_b=vector(0,ncols-1);
	double *vhat_temp;
	vhat_temp =vector(0,ncols-1);

	for (int iteration=0;iteration<max_iter;iteration++)  
	{	//main iteration loop
		
		



		//bit-to-check messages
		for ( i=0;i<mrows;i++)
		{
			for (j=0;j<max_check_degree;j++)
			{

				col_index=check_node_ones_matrix[i][j];

				if (col_index==BIGVALUE_COLS)					// if element is value put in to fill zeros, then break
				{
					break;
				}

				bitmessage_temp[i][(int)col_index]=tanh (   (     -(sg_array[i][(int)col_index]) + sa_array[i][(int)col_index]     ) / 2  );
			}

		}


		//check-to-bit-messages
		for (int u=0;u<mrows;u++)
		{//across mrows
			double temp=1;
			
			for (int v =0;v<max_check_degree;v++)	
			{//across columns
				
				element=check_node_ones_matrix[u][v];  // this is value of column in H which has '1' for this row

				if (element==BIGVALUE_COLS)					// if element is value put in to fill zeros, then break
				{
					break;
				}
				
				for (int w=0;w<max_check_degree;w++)
					{//accross columns again
						
						col_index=check_node_ones_matrix[u][w];

						if(check_node_ones_matrix[u][w]!=element && check_node_ones_matrix[u][w]!=BIGVALUE_COLS)
						{//second if
							temp=temp* bitmessage_temp[u][(int)col_index];
						}//second if

					}//accross columns again

					sa_array[u][(int)element]=-2*atanh(temp);
					temp=1;
				
			
			}//across columns
		
		}//accross mrows    



		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		// to print to screen the 'sa_array' matrix
	/*	for ( i=0;i<mrows;i++)
		{
			mexPrintf("\n");
			for(int j=0;j<ncols;j++)
			{
				mexPrintf("\t%f",sa_array[i][j]);
			}
			mexPrintf("\n");
		} */
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




		//sum across columns
		for ( i=0;i<ncols;i++)
		{
			double temp=0;
			for(int j=0;j<max_variable_degree;j++)
			{	
				row_index=variable_node_ones_matrix[j][i];
				if(row_index==BIGVALUE_ROWS)
					break;


				temp=temp+sa_array[(int)row_index][i];//sa_array is equal to old bitmessage_2
			}
		
			sum_of_b[i]=temp;

		}


	
		//update the APP LLR and hard decision
		for ( i=0;i<ncols;i++)
		{
			for(int j=0;j<max_variable_degree;j++)
			{
				
				row_index=variable_node_ones_matrix[j][i];
				if(row_index==BIGVALUE_ROWS)
					break;

				sg_array[(int)row_index][i]=sum_of_b[i]+gamma_n[i];
		
	
					//mexPrintf("\t%f",sg_array[j][i]);
				
			}

			vhat_temp[i]=sum_of_b[i]+gamma_n[i];
		
			//hard decision
			if (vhat_temp[i]>=0)
				*(vhat+i)=1;
			else
				*(vhat+i)=0;
			//mexPrintf("\n\t vhat is :%f\t",*(vhat+i));
		
		}
		


		//check if valid codeword
	
	    int parity=0,cumsum=0;

		for ( j=0;j<mrows;j++)
		{
			
			for (i=0;i<ncols;i++)	//multiply codeword with row of H
			{
				*(vhat_temp+i)=H[j][i]* *(vhat+i);
			}
			
			parity=0;
			for(i=0;i<ncols;i++)
			{
				parity=(parity + (int)(*(vhat_temp+i)) )%2;
				 
			}
		
			//mexPrintf("\n\tparity is :%d\t",parity);
			cumsum=parity+cumsum;
			//mexPrintf("\n\tInput cumsum is :%f\t",cumsum);
			
			if (cumsum==1) //will happen when parity is 1: means not valid codeword, so continue with iterations
			{	break;	}
		}

		if (cumsum==0) //valid codeword
		{	
			*iter_out=iteration+1;
			return; 
		}

	}	//main iteration loop
		
	*iter_out=iteration;

}//function braces



double atanh(double x) {
  double epsilon;	
  epsilon = pow(10,-16);
  
  if(x>(1-epsilon)) return INF;
  if(x<(-1+epsilon)) return -INF;
  return 0.5*log((1+x)/(1-x));
}



//      new              0        1          2              3                 4                5                  6                 
//vhat=decode_ldpc_new(max_iter,gamma_n,check_node_ones,max_check_degree,BIGVALUE_COLS-1,variable_node_ones,max_variable_degree,

//     7            8   9
//BIGVALUE_ROWS-1,rows,cols);



void mexFunction( int nlhs, mxArray *plhs[], 
				  int nrhs, const mxArray*prhs[] )
{
	double *vhat, *iter_out,*gamma_n, *check_node_ones, *variable_node_ones; /*pointer variables for input Matrices*/
	double max_iter,max_check_degree,max_variable_degree, BIGVALUE_COLS,BIGVALUE_ROWS; 
	int mrows,ncols;


	gamma_n  = mxGetPr(prhs[1]); //pointer to initial APP LLR
	check_node_ones=mxGetPr(prhs[2]);  //pointer to matrix containing column indeces in each row which are '1'
	max_check_degree=mxGetScalar(prhs[3]); //the maximum check degree
	BIGVALUE_COLS=mxGetScalar(prhs[4]);

	variable_node_ones=mxGetPr(prhs[5]); //pointer to matrix containing row indeces in each column which are '1'
	max_variable_degree=mxGetScalar(prhs[6]); //the maximum variable degree
	BIGVALUE_ROWS=mxGetScalar(prhs[7]);



	//No = mxGetScalar(prhs[0]);       // value of No
	max_iter = mxGetScalar(prhs[0]); // maximum iterations
	
//	mexPrintf("\n\tInput Arg No is :%f\t",No);
//mexPrintf("\n\tInput Arg max_iter is :%f\t",max_iter);

	  mrows = mxGetScalar(prhs[8]); // number of rows of H)
	  ncols = mxGetScalar(prhs[9]); // number of cols of H)

	 plhs[0] = mxCreateDoubleMatrix(1,ncols, mxREAL); /*matrix for output*/
	 vhat = mxGetPr(plhs[0]);	/*pointer to output*/
	 plhs[1] = mxCreateDoubleScalar(0);
	 iter_out = mxGetPr(plhs[1]);

	 decode(max_iter,vhat,mrows,ncols,iter_out,gamma_n,check_node_ones,max_check_degree,BIGVALUE_COLS,variable_node_ones,max_variable_degree,BIGVALUE_ROWS);


}