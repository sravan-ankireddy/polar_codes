#include	<stdio.h>
#include	<stdlib.h>
#include	<math.h>
#include	<time.h>



int *ivector(int ilow,int ihigh);
int **imatrix(int ilow,int ihigh,int jlow,int jhigh);
double *vector(int ilow,int ihigh);
double **matrix(int ilow,int ihigh,int jlow,int jhigh);



double *vector(int ilow,int ihigh)
/*************************************************************************
*   This routine allocates space for a double vector using malloc and sets 
*   the offset so that the array can be referenced from ilow to ihigh.
*
*   ALGORITHM:  taken from Numerical Recipes in C,Press et.al.,
*               Cambridge Press 1990,pg705.
*   INPUTS: ilow - lowest value in index range
*           ihigh - highest value in index range
*   OUTPUTS: none
*   RETURNS: pointer to allocated space with the proper offset 
*
*   K. Chugg 10/15/90 
**************************************************************************/
{
    double *v;
	v=(double *)mxCalloc((size_t) (ihigh-ilow+1), sizeof(double));  // to initialize elements to '0'
    //if (!v) nrerror("allocation failure in vector()");
    return v-ilow;
}


double **matrix(int ilow,int ihigh,int jlow,int jhigh)
/*************************************************************************
*   This routine allocates space for a double matrix using malloc and sets 
*   the offset so that the matrix can be referenced from ilow to ihigh and
*   jlow to jhigh.
*
*   ALGORITHM:  taken from Numerical Recipes in C,Press et.al.,
*               Cambridge Press 1990,pg706.
*   INPUTS: ilow - lowest value in row range
*           ihigh - highest value in row range
*           jlow - lowest value in col range
*           jhigh - highest value in col range
*           
*   OUTPUTS: none
*   RETURNS: pointer to pointers to allocated space with the proper offsets 
*   NOTES: unlike with a 2-D array m[i] represents a "vector" - the ith row
*   K. Chugg 10/21/90 
**************************************************************************/
{
    int i;
    double **m;
    
    /* Allocate pointers to the rows  */
	m=(double **)mxCalloc((size_t) (ihigh-ilow+1), sizeof(double*));  // to initialize elements to '0'
    //if (!m) nrerror("allocation error in matrix()");
    m -= ilow;
    
    /* allocate space for the rows and set the pointers to them  */
    for (i=ilow;i<=ihigh;i++)
        {
		m[i]=(double *)mxCalloc((size_t) (jhigh-jlow+1), sizeof(double));  // to initialize elements to '0'
       // if(!m[i]) nrerror("allocation error 2 in dmatrix()");
        m[i] -= jlow;
        }
    return m;
}


int *ivector(int ilow,int ihigh)
/*************************************************************************
*   This routine allocates space for a int vector using malloc and sets the
*   offset so that 
*   the array can be referenced from ilow to ihigh.
*
*   ALGORITHM:  taken from Numerical Recipes in C,Press et.al.,
*               Cambridge Press 1990,pg705.
*   INPUTS: ilow - lowest value in index range
*           ihigh - highest value in index range
*   OUTPUTS: none
*   RETURNS: pointer to allocated space with the proper offset 
*
*   K. Chugg 10/15/90 
**************************************************************************/
{
    int *v;

    v=(int *)mxMalloc((size_t) ((ihigh-ilow+1)*sizeof(int)));
    //if (!v) nrerror("allocation failure in vector()");
    return v-ilow;
}

int **imatrix(int ilow,int ihigh,int jlow,int jhigh)
/*************************************************************************
*   This routine allocates space for a int matrix using malloc and sets the
*   offset so that the matrix can be referenced from ilow to ihigh and 
*   jlow to jhigh.
*
*   ALGORITHM:  taken from Numerical Recipes in C,Press et.al.,
*               Cambridge Press 1990,pg706.
*   INPUTS: ilow - lowest value in row range
*           ihigh - highest value in row range
*           jlow - lowest value in col range
*           jhigh - highest value in col range
*           
*   OUTPUTS: none
*   RETURNS: pointer to pointers to allocated space with the proper offsets 
*   NOTES: unlike with a 2-D array m[i] represents a "vector" - the ith row
*   K. Chugg 10/21/90 
*************************************************************************/
{
    int i;
    int **m;
    
    /* Allocate pointers to the rows  */
    m=(int **)mxMalloc((size_t) (ihigh-ilow+1)*sizeof(int*));
    //if (!m) nrerror("allocation error in matrix()");
    m -= ilow;
    
    /* allocate space for the rows and set the pointers to them  */
    for (i=ilow;i<=ihigh;i++)
        {
        m[i]=(int *)mxMalloc((size_t) (jhigh-jlow+1)*sizeof(int));
        //if(!m[i]) nrerror("allocation error 2 in dmatrix()");
        m[i] -= jlow;
        }
    return m;
}