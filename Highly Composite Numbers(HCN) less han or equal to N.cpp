
/// Highly Composite Numbers are product of Primorials (Primorials are same as Factorials, but instead of natural number we multiply primes.
/// p3=2×3×5 and p5=2×3×5×7×11 )

#include<bits/stdc++.h>
using namespace std;

int prime[] = {2, 3, 5, 7, 11, 13, 17, 19, 23 };
/// It contains first 9 primes so it will work correctly for N≤109

int resNum, resDiv, n ;
void recur ( int pos, int limit, long long num, int div )
{
    if ( div > resDiv )   /// Get the number with highest NOD
    {
        resNum = num;
        resDiv = div;
    }
    else if ( div == resDiv && num < resNum )   /// In case of tie, take smaller number
    {
        resNum = num;
    }

    if ( pos == 9 )
        return; /// End of prime list

    long long p = prime[pos];

    for ( int i = 1; i <= limit; i++ )
    {
        if ( num * p > n )
            break;
        recur ( pos + 1, i, num * p, div * ( i + 1 ) );
        p *= prime[pos];
    }
}
int main()
{
    n = 1000000000;

    resNum = resDiv = 0 ;
    recur(0,30,1,1);

    /// limit parameter is set to 30 since 2^30 > N .

    printf("%d %d\n",resNum,resDiv);
}


/// For programming contest, we could memorize values of HCN that comes frequently.
/// Mainly 1344 for N≤10^9 and 103,680 for N≤10^18.




unsigned long long n, res;
int p, primes[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 51, 53, 59, 61, 67, 71};

unsigned long long mul(unsigned long long a, unsigned long long b){
    unsigned long long res = 0;
    while (b){
        if (b & 1LL) res = (res + a);
        if (res >= n) return 0;
        a = (a << 1LL);
        b >>= 1LL;
    }
    return res;
}

void backtrack(int i, int lim, unsigned long long val, unsigned long long r){
    if (r > res) res = r;
    if (i == p) return;

    int d;
    unsigned long long x = val;

    for (d = 1; d <= lim; d++){
        x = mul(x, primes[i]);
        if (x == 0) return;
        backtrack(i + 1, d, x, r * (d + 1));
    }
}

int main(){
    /* Tested for n <= 10^18 */

    p = sizeof(primes) / sizeof(int);

    while (scanf("%llu", &n) != EOF){
        res = 0;
        backtrack(0, 100, 1, 1);
        printf("Maximum number of divisors of any number less than %llu = %llu\n", n, res);
    }
    return 0;
}

