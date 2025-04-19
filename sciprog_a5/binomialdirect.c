#include <stdio.h>

// the intermediate factorials can become very large very fast
// for n>=12, the function overflows, i.e. starts to give garbage values
// overflow happens because for the datea type int C uses 4 bytes,
// which translates to the value range of -2,147,483,648 to 2,147,483,647.
// so instead we could use data types with bigger sizes like long long or unsigned long long

// from the worksheet
int factorial(int n)
{
    if (n > 1)
    {
        return n * factorial(n - 1);
    }
    else
    {
        return 1;
    }
}

int binomialdirect(int n, int k)
{
    return factorial(n) / (factorial(k) * factorial(n - k));
}

int main()
{
    int n, k;

    printf("enter n and k: ");
    scanf("%d %d", &n, &k);

    int result = binomialdirect(n, k);

    printf("binomial coefficient (n choose k) is: %d\n", result);

    return 0;
}
