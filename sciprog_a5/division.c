#include <stdio.h>

// recursive integer division without remainder
int division(int m, int n)
{
    if (m < n)
    {
        return 0;
    }
    else
    {
        return 1 + division(m - n, n);
    }
}

int main()
{
    int m, n;

    printf("enter m and n (m >= 0, n > 0): ");
    scanf("%d %d", &m, &n);

    // basic check
    if (n <= 0)
    {
        printf("invalid input: n must be > 0\n");
        return 1;
    }

    int result = division(m, n);
    printf("result of m / n (integer division) is: %d\n", result);

    return 0;
}
