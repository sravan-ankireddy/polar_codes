//  Copyright (c) 2005 by Shaikh Faisal Zaheer, faisal.zaheer@gmail.com
//  $Revision: 1.0 $  $Date: 10/11/2005 $


/*LDPC Decoder*/

#include "mex.h"
#include "matrix.h"			// for Matlab mx and mex fuctions
#include "math.h"
#include <stdlib.h>			// what for 
#include "decodeutil.h"

#define INF 50  // maximum value of extrinsic LLR checknode to bitnode messages

double atanh(double x);

void decode(double No,double *h,double max_iter,
			double *sg, double *vhat, int mrows, int ncols, double *iter_out, double *gamma_n )

{//function braces

	
	double **sg_array,**sa_array, **bitmessage_temp ;
	double **H;
	
	sg_array		= matrix(0,mrows-1,0,ncols-1);	
	sa_array		= matrix(0,mrows-1,0,ncols-1);	
	H				= matrix(0,mrows-1,0,ncols-1);	

		int i=0, j=0;

		for ( i=0;i<ncols;i++)
		{
			for(j=0;j<mrows;j++)
			{
				sg_array[j][i]=*(sg++);
				sa_array[j][i]=0;
				H[j][i]=*(h++);
			}
		}

	/*	for ( i=0;i<mrows;i++)
		{
			mexPrintf("\n");
			for(int j=0;j<ncols;j++)
			{
				mexPrintf("\t%f",H[i][j]);
			}
			mexPrintf("\n");
		} */

	bitmessage_temp = matrix(0,mrows-1,0,ncols-1);
	double *sum_of_b;
	sum_of_b=vector(0,ncols-1);
	double *vhat_temp;
	vhat_temp =vector(0,ncols-1);

	for (int iteration=0;iteration<max_iter;iteration++)  
	{	//main iteration loop
		
		//bit-to-check messages
		for ( i=0;i<ncols;i++)
		{
			for(int j=0;j<mrows;j++)
			{
				if (H[j][i]==1)
				bitmessage_temp[j][i]=tanh (   (     -(sg_array[j][i]) + sa_array[j][i]     ) / 2  );
			}
		}

		//check-to-bit-messages
		for (int u=0;u<mrows;u++)
		{//across mrows
			double temp=1;
			
			for (int v =0;v<ncols;v++)	
			{//across columns
				
				if (H[u][v]==1)
				{// first if
					
					for (int w=0;w<ncols;w++)
					{//accross columns again
						
						if(H[u][w]==1 && v!=w)
						{//second if
							temp=temp* bitmessage_temp[u][w];
						}//second if

					}//accross columns again

					sa_array[u][v]=-2*atanh(temp);
					temp=1;
				}//first if
			
			}//across columns
		
		}//accross mrows


		//sum across columns
		for ( i=0;i<ncols;i++)
		{
			double temp=0;
			for(int j=0;j<mrows;j++)
			{
				temp=temp+sa_array[j][i];//sa_array is equal to old bitmessage_2
			}
		
			sum_of_b[i]=temp;
		}

		//
		for ( i=0;i<ncols;i++)
		{
			for(int j=0;j<mrows;j++)
			{
				if (H[j][i]==1)
				{
					sg_array[j][i]=sum_of_b[i]+gamma_n[i];
					vhat_temp[i]=sum_of_b[i]+gamma_n[i];
				
					//mexPrintf("\t%f",sg_array[j][i]);
				}
			}
		}




		//hard decision
		for (i=0;i<ncols;i++)
		{
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






	//								0 1    2      3     4          
//[vhat,iteration]=decode_ldpc_log(No,h,max_iter,sg,gamma_n)
void mexFunction( int nlhs, mxArray *plhs[], 
				  int nrhs, const mxArray*prhs[] )
{
	double *h,*vhat,*sg, *iter_out,*gamma_n; /*pointer variables for input Matrices*/
	double No,max_iter; 
	int mrows,ncols;


	h		= mxGetPr(prhs[1]); /*pointer to H */
	sg		= mxGetPr(prhs[3]); //pointer to sg
	gamma_n  = mxGetPr(prhs[4]); //pointer to initial APP LLR


	No = mxGetScalar(prhs[0]);       // value of No
	max_iter = mxGetScalar(prhs[2]); // maximum iterations
	
//	mexPrintf("\n\tInput Arg No is :%f\t",No);
//mexPrintf("\n\tInput Arg max_iter is :%f\t",max_iter);

	  mrows = mxGetM(prhs[1]); // number of rows of H)
	  ncols = mxGetN(prhs[1]); // number of cols of H)

	 plhs[0] = mxCreateDoubleMatrix(1,ncols, mxREAL); /*matrix for output*/
	 vhat = mxGetPr(plhs[0]);	/*pointer to output*/
	 plhs[1] = mxCreateDoubleScalar(0);
	 iter_out = mxGetPr(plhs[1]);

	 decode(No,h,max_iter,sg,vhat,mrows,ncols,iter_out,gamma_n);


}