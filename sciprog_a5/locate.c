#include <stdio.h>

int main()
{
    double L;
    double x;
    double y;

    printf("Write the number L for the square: ");
    scanf("%lf", &L);

    printf("Write the x coordinate of the point: ");
    scanf("%lf", &x);

    printf("Write the y coordinate of the point: ");
    scanf("%lf", &y);

    if (x > 0 && x < L && y > 0 && y < L)
    {
        printf("The point (%.2lf, %.2lf) is inside the square given by (0,0), (%.2lf,0), (%.2lf,%.2lf), (0,%.2lf)\n", x, y, L, L, L);
    }
    else if ((x == 0 || x == L) && (y >= 0 && y <= L) || (y == 0 || y == L) && (x >= 0 && x <= L))
    {
        printf("The point (%.2lf, %.2lf) is on the boundary of the square.\n", x, y);
    }
    else
    {
        printf("The point (%.2lf, %.2lf) is outside the square given by (0,0), (%.2lf,0), (%.2lf,%.2lf), (0,%.2lf)\n", x, y, L, L, L);
    }

    return 0;
}
