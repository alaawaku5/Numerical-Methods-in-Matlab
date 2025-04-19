#include <stdio.h>

// checking collinearity using determinant: 1/2
void points(double x, double y, double u, double v, double a, double b)
{
    double D = x * (v - b) + u * (b - y) + a * (y - v);

    if (D == 0)
    {
        printf("The points lie on the same line.\n");
    }
    else
    {
        printf("The points do NOT lie on the same line.\n");
    }
}

int main()
{
    double x, y, u, v, a, b;

    printf("Enter the first point (x y): ");
    scanf("%lf %lf", &x, &y);

    printf("Enter the second point (u v): ");
    scanf("%lf %lf", &u, &v);

    printf("Enter the third point (a b): ");
    scanf("%lf %lf", &a, &b);

    points(x, y, u, v, a, b);

    return 0;
}
