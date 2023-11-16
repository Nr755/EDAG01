#include <stdio.h>

int s;
int a;

int main()
{
    s = "";
    a = 3;
    while (a != 1)
    {
        s += toString(a);
        if (a % 2 == 0)
            a = a / 2;
        else
            a = 3 * a + 1;
    }

    print(s);
}