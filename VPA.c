#include <stdio.h>

int main(){

    printf("Size of FLOAT type (bytes): %d\n", sizeof(float));
    printf("EPSILON value (FLOAT): %.5e\n", __FLT_EPSILON__);
    printf("\n");

    printf("Size of DOUBLE type (bytes): %d\n", sizeof(double));
    printf("EPSILON value (FLOAT): %.5le\n", __DBL_EPSILON__);
    printf("\n");

    printf("Size of LONG DOUBLE type (bytes): %d\n", sizeof(long double));
    printf("EPSILON value (FLOAT): %.5Le\n", __LDBL_EPSILON__);
    printf("\n");

    return 0;
}