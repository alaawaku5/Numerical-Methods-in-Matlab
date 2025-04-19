#include <stdio.h>

int main()
{
    double x;
    int n;

    printf("amount of money won (x>0):");
    scanf("%lf", &x);
    // we expect to receive a "%lf" -> a double
    // and want to put it at x's address, where & specifies address,

    printf("number of friends (n > 0): ");
    scanf("%d", &n);

    double y = x / (n + 1);
    printf("each of you will get %.2lf euros.\n", y);

    return 0;
}