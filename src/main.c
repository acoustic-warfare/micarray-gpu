#include "beam.h"

// Dummy param for creating a simple array of BYTES size
#define BYTES 64

extern void cudaBeamWrapper(const int *res, const int *first, const int *last, int n_bytes);

// int main(int argc, char *argv[])

// {

//     int *a = (int *)malloc(sizeof(int) * 2);

//     a[0] = 2;

//     a[1] = 3;

//     printf("a[0]: %d, a[1]: %d\n", a[0], a[1]);

//     kernel_wrapper(a);

//     printf("a[0]: %d, a[1]: %d\n", a[0], a[1]);

//     free(a);

//     return 0;
// }

int main(int argc, char **argv)
{
    const int arraySize = BYTES;
    // const int res[arraySize] = {0};
    int res[arraySize];
    memset(res, 0, arraySize * sizeof(int));
    int first[arraySize];
    int last[arraySize];

    // Inititate random values
    int i;
    for (i = 0; i < BYTES; i++)
    {
        // first[i] = rand();
        first[i] = i;
    }
    for (i = 0; i < BYTES; i++)
    {
        // last[i] = rand();
        last[i] = 2;
    }

    cudaBeamWrapper(res, first, last, arraySize);
    int loop;
    for (loop = 0; loop < BYTES; loop++)
        printf("%d ", res[loop]);
    printf("\n");

    return 0;
}