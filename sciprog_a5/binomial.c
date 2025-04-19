#include <stdio.h>

int binomial(int n, int k)
{
    if (k == 0 || k == n)
    {
        return 1;
    }
    else
    {
        return binomial(n - 1, k - 1) + binomial(n - 1, k);
    }
}

int main()
{
    int n, k;

    // ask user for input
    printf("enter n and k: ");
    scanf("%d %d", &n, &k);

    // basic check (optional but helpful)
    if (k < 0 || k > n)
    {
        printf("invalid input: make sure 0 <= k <= n\n");
        return 1;
    }

    int result = binomial(n, k);
    printf("binomial coefficient (n choose k) is: %d\n", result);

    return 0;
}
