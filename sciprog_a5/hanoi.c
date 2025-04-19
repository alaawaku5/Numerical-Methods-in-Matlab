#include <stdio.h>

// Recursive function to solve Tower of Hanoi
void hanoi(int m, int i, int j)
{
    if (m == 1)
    {
        printf("Move a disk from rod %d to rod %d.\n", i, j);
    }
    else
    {
        int k = 6 - i - j;                                    // the third rod
        hanoi(m - 1, i, k);                                   // Step 1: move m-1 disks from i to k
        printf("Move a disk from rod %d to rod %d.\n", i, j); // Step 2
        hanoi(m - 1, k, j);                                   // Step 3: move m-1 disks from k to j
    }
}

int main()
{
    int n;
    printf("Enter the number of disks: ");
    scanf("%d", &n);

    hanoi(n, 1, 3); // Move n disks from rod 1 to rod 3

    return 0;
}

// # hanoi(5, 1, 3)
// ├── hanoi(4, 1, 2)
// │   ├── hanoi(3, 1, 3)
// │   │   ├── hanoi(2, 1, 2)
// │   │   │   ├── hanoi(1, 1, 3)
// │   │   │   └── hanoi(1, 3, 2)
// │   │   └── hanoi(2, 2, 3)
// │   │       ├── hanoi(1, 2, 1)
// │   │       └── hanoi(1, 1, 3)
// │   └── hanoi(3, 3, 2)
// │       ├── hanoi(2, 3, 1)
// │       │   ├── hanoi(1, 3, 2)
// │       │   └── hanoi(1, 2, 1)
// │       └── hanoi(2, 1, 2)
// │           ├── hanoi(1, 1, 3)
// │           └── hanoi(1, 3, 2)
// ├── PRINT (move largest disk 5 from 1 to 3)
// └── hanoi(4, 2, 3)
//     ├── hanoi(3, 2, 1)
//     │   ├── hanoi(2, 2, 3)
//     │   │   ├── hanoi(1, 2, 1)
//     │   │   └── hanoi(1, 1, 3)
//     │   └── hanoi(2, 3, 1)
//     │       ├── hanoi(1, 3, 2)
//     │       └── hanoi(1, 2, 1)
//     └── hanoi(3, 1, 3)
//         ├── hanoi(2, 1, 2)
//         │   ├── hanoi(1, 1, 3)
//         │   └── hanoi(1, 3, 2)
//         └── hanoi(2, 2, 3)
//             ├── hanoi(1, 2, 1)
//             └── hanoi(1, 1, 3)
