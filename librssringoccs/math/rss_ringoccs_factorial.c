/*  Header file which contains aliases for the function in the standard C     *
 *  library math.h. This allows compatibility of C89 and C99 math.h headers.  */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

long rssringoccs_Factorial(int N)
{
    int n;
    long factorial = 1;

    for (n=1; n<=N; ++n)
        factorial *= n;

    return factorial;
}

long rssringoccs_Falling_Factorial(int x, int N)
{
    int n;
    long falling_factorial = x;

    if (N == 0)
        return 1;
    else if (N > x)
        return 0;
    else
    {
        for (n=1; n<N; ++n)
            falling_factorial *= (x-n);
        return falling_factorial;
    }
}
