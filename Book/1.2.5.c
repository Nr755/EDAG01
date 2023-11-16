#include <stdbool.h>
#include <stdio.h>

bool leap_year(unsigned short int year)
{
    if (year % 4 != 0)
    {
        return false;
    }
    else if (year % 100 != 0)
    {
        return true;
    }
    else if (year % 400 != 0)
    {
        return false;
    }
    else
    {
        return true;
    }
}

int main()
{
    unsigned short int year;
    printf("Enter a year: ");
    scanf("%hu", &year);
    if (leap_year(year))
    {
        printf("%hu is a leap year\n", year);
    }
    else
    {
        printf("%hu is not a leap year\n", year);
    }
    return 0;
}
