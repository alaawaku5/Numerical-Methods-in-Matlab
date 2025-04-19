#include <stdio.h>
#include <math.h>

int main()
{
    double x, y, z;
    double max_val, min_val, mid_val;

    printf("Enter three numbers (x, y, z):\n");
    scanf("%lf %lf %lf", &x, &y, &z);

    max_val = fmax(x, fmax(y, z));
    min_val = fmin(x, fmin(y, z));

    mid_val = (x + y + z) - max_val - min_val;

    // Output
    printf("Numbers in descending order: %.2lf %.2lf %.2lf\n", max_val, mid_val, min_val);

    return 0;
}
