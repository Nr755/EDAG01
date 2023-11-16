#include <stdio.h>
#include <math.h>

int main()
{
    double x[] = {0, 0.1, 0.2, 0.3};
    int n = sizeof(x) / sizeof(x[0]);

    printf("|   x   | sqrt(x) |  e^x   |\n");
    printf("|-------|---------|--------|\n");

    for (int i = 0; i < n; i++)
    {
        double sqrt_x = sqrt(x[i]);
        double exp_x = exp(x[i]);
        printf("| %5.1f | %7.4f | %6.4f |\n", x[i], sqrt_x, exp_x);
    }

    return 0;
}
