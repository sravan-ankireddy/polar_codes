#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include <stdint.h>
#include <inttypes.h>
#include <arm_neon.h>


int main()
{


	int8x16_t a, b;

	uint8x16_t result;

	int8_t a_int[16] = {1,2,-3,4,5,-6,7,8,-1,2,3,4,-5,6,7,-8};
	int8_t b_int[16] = {5,6,-7,3,4,-9,22,3,-6,7,78,1,2,3,4,9};
	uint8_t result_int[16];
	uint8_t result_neon[16];

	int i;

	for(i = 0; i < 16; i++)
	{
		if (a_int[i] >= b_int[i])
			result_int[i] = 0;
		else
			result_int[i] = 1;

		printf("%d\t", result_int[i]);
	}

	printf("\n");

	a = vld1q_s8(a_int);
	b = vld1q_s8(b_int);

	result = vcgeq_s8(b,a);

	vst1q_u8(result_neon,result);

	for(i = 0; i < 16; i++)
	{
		printf("%d\t", result[i]);
	}
	printf("\n");
	for(i = 0; i < 16; i++)
	{
		printf("%d\t", result_neon[i]);
	}

	printf("\nxor test\n");
	printf("%d\t", 255^0);
	printf("%d\n", 1^0);

	printf("%d\t", 255^255);
	printf("%d\t", 1^1);
	printf("%d\t", 1^255);
	printf("%d\n", 255^1);
}