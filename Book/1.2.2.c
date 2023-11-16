
#include <stdio.h>

int main()
{
    unsigned short us = ~0;
    unsigned int ui = ~0;
    unsigned long ul = ~0;
    unsigned long long ull = ~0;

    int bits_us = 0;
    while (us)
    {
        bits_us++;
        us >>= 1;
    }

    int bits_ui = 0;
    while (ui)
    {
        bits_ui++;
        ui >>= 1;
    }

    int bits_ul = 0;
    while (ul)
    {
        bits_ul++;
        ul >>= 1;
    }

    int bits_ull = 0;
    while (ull)
    {
        bits_ull++;
        ull >>= 1;
    }

    printf("Number of bits for unsigned short: %d\n", bits_us);
    printf("Number of bits for unsigned int: %d\n", bits_ui);
    printf("Number of bits for unsigned long: %d\n", bits_ul);
    printf("Number of bits for unsigned long long: %d\n", bits_ull);

    return 0;
}
