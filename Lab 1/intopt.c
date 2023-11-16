#include <stdio.h>
#include <stdlib.h>

int main()
{
    int m, n;
    scanf("%d %d", &m, &n);

    double **matrix = (double **)malloc(m * sizeof(double *));
    double *coefficients = (double *)malloc(m * sizeof(double));
    double *values = (double *)malloc(m * sizeof(double));

    int i, j;

    // Read coefficients
    for (i = 0; i < m; i++)
    {
        scanf("%lf", &coefficients[i]);
    }

    // Read matrix elements
    for (i = 0; i < m; i++)
    {
        matrix[i] = (double *)malloc(n * sizeof(double));
        for (j = 0; j < n; j++)
        {
            scanf("%lf", &matrix[i][j]);
        }
    }

    // Read values
    for (i = 0; i < m; i++)
    {
        scanf("%lf", &values[i]);
    }

    printf("max Z = ");
    for (i = 0; i < n; i += 1)
    {
        if (i > 0)
        {
            printf("+ %10.3lfx_%d ", coefficients[i], i);
        }
        else
        {
            printf("%10.3lfx_%d ", coefficients[i], i);
        }
    }
    printf("\n");

    for (i = 0; i < m; i += 1)
    {
        printf("\t");
        for (j = 0; j < n; j += 1)
        {
            if (j > 0)
            {
                printf("+ %10.3lfx_%d ", matrix[i][j], j);
            }
            else
            {
                printf("%10.3lfx_%d ", matrix[i][j], j);
            }
        }
        printf(" \u2264 %10.3lf\n", values[i]);
        // printf("\n");
    }

    // Free dynamically allocated memory
    free(values);
    free(coefficients);
    // for (i = 0; i < m; i++)
    // {
    //     free(matrix[i]);
    // }
    free(matrix);

    return 0;
}
