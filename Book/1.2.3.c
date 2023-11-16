
#include <stdio.h>
#include <time.h>

int main()
{
    clock_t start_time = clock();
    double sum = 0.0;
    int terms = 0;
    while (sum <= 18.0)
    {
        terms++;
        sum += 1.0 / terms;
    }
    clock_t end_time = clock();
    double time_used = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;
    printf("Number of terms: %d\n", terms);
    printf("Sum: %.15lf\n", sum);
    printf("Time used: %.15lf seconds\n", time_used);
    return 0;
}
