#include <stdio.h>

double scalarproduct(double a, double b, double c, double x, double y, double z)
{
    return a * x + b * y + c * z;
}

int main()
{
    double a, b, c; // first vector
    double x, y, z; // second vector

    printf("enter a, b, c (first vector): ");
    scanf("%lf %lf %lf", &a, &b, &c);

    printf("enter x, y, z (second vector): ");
    scanf("%lf %lf %lf", &x, &y, &z);

    double result = scalarproduct(a, b, c, x, y, z);

    printf("the scalar product is: %.2lf\n", result);

    return 0;
}
