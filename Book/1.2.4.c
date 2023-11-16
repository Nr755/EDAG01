
#include <stdio.h>

int hex_digit_sum(unsigned int num)
{
    int sum = 0;
    while (num != 0)
    {
        sum += num & 0xf;
        num >>= 4;
    }
    return sum;
}

int main()
{
    unsigned int num = 0x12345678;
    int sum = hex_digit_sum(num);
    printf("Hexadecimal digit sum of %x is %d\n", num, sum);
    return 0;
}
