#include <stdio.h>

// recursive fibonacci function
int fibonacci(int n)
{
    if (n == 0)
    {
        return 0;
    }
    if (n == 1)
    {
        return 1;
    }
    else
    {
        return fibonacci(n - 1) + fibonacci(n - 2);
    }
}

// fibonacci(5)
// ├── fibonacci(4)
// │   ├── fibonacci(3)
// │   │   ├── fibonacci(2)
// │   │   │   ├── fibonacci(1)
// │   │   │   └── fibonacci(0)
// │   │   └── fibonacci(1)
// │   └── fibonacci(2)
// │       ├── fibonacci(1)
// │       └── fibonacci(0)
// └── fibonacci(3)
//     ├── fibonacci(2)
//     │   ├── fibonacci(1)
//     │   └── fibonacci(0)
//     └── fibonacci(1)

int main()
{
    int n;

    printf("enter n (non-negative): ");
    scanf("%d", &n);

    if (n < 0)
    {
        printf("invalid input: n must be >= 0\n");
        return 1;
    }

    int result = fibonacci(n);
    printf("fibonacci number at position %d is: %d\n", n, result);

    return 0;
}

// int fibonacci(int n)
// {
//     if (n == 0) return 0;
//     if (n == 1) return 1;

//     int a = 0, b = 1, temp;
//     for (int i = 2; i <= n; i++)
//     {
//         temp = a + b;
//         a = b;
//         b = temp;
//     }
//     return b;
// }

// vars used:
// - a: stores F(n-2)
// - b: stores F(n-1)
// - temp: stores F(n) = a + b

// at each step i from 2 to n:
//     temp = a + b     // compute the next Fibonacci number
//     a = b            // shift a to F(n-1)
//     b = temp         // update b to the new F(n)

// example for n = 5:
// Initial: a = 0 (F(0)), b = 1 (F(1))
// Step 2: temp = 0 + 1 = 1 → a = 1, b = 1
// Step 3: temp = 1 + 1 = 2 → a = 1, b = 2
// Step 4: temp = 1 + 2 = 3 → a = 2, b = 3
// Step 5: temp = 2 + 3 = 5 → a = 3, b = 5

// returns b = 5, which is F(5)

// runtime: O(n)
// space complexity: O(1)