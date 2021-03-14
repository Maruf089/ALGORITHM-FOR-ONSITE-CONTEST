/**
    Dependencies:
        1. MillerRabin
    How to Use it?
        1. Call pollardRho.clear();
        2. Call pollardRho.getPrimeFactorization(n);
    See sample main() function below
*/
#include<bits/stdc++.h>
using namespace std;
#define ll long long int


typedef long long vlong;

class MillerRabin {
    private:

    /** https://miller-rabin.appspot.com/ confirms that the following base covers 2^64**/
    vlong prime[7] = { 2, 325, 9375, 28178, 450775, 9780504, 1795265022 };
    int psize = 7;

    vlong bigmod ( __int128 a, __int128 p, vlong mod ) {
        __int128 x = a % mod, res = 1;
        while ( p ) {
            if ( p & 1 ) res = res * x % mod;
            x = x * x % mod;
            p >>= 1;
        }
        return res;
    }

    ///Witness to compositeness of n
    ///n - 1 = 2^s * d, where d is odd
    bool witness ( vlong a, vlong d, vlong s, vlong n ) {
        __int128 r = bigmod( a, d, n );
        if ( r == 1 || r == n - 1 ) return false;
        int i;
        for ( i = 0; i < s - 1; i++ ) {
            r = r * r % n;
            if ( r == 1 ) return true;
            if ( r == n - 1 ) return false;
        }
        return true;
    }

    public:
    bool isPrime ( vlong n ) {
        if ( n <= 1 ) return false;

        vlong p = n - 1, s = 0;
        while ( ! ( p & 1 ) ) {
            p /= 2;
            s++;
        }
        vlong d = p;
        p = n - 1;

        for ( int i = 0; i < psize && prime[i] < n; i++ ) {
            if ( witness ( prime[i], d, s, n ) ) return false;
        }
        return true;
    }
} millerRabin;


class PollardRho {
    private:

    MillerRabin millerRabin;

    int prime[50000], status[50000], primeSize;
    void sieve() {
        primeSize = 0;
        memset( status, 0, sizeof status );

        status[0] = status[1] = 1;
        int n = 46340;
        for ( int i = 4; i <= n; i += 2 ) status[i] = 1;

        int sqrtn = sqrt(n);
        for ( int i = 3; i <= sqrtn; i += 2 ){
            for ( int j = i*i; j <= n; j += 2 * i ) {
                status[j] = 1;
            }
        }

        prime[primeSize++] = 2;
        for ( int i = 3; i <= n; i += 2 ) {
            if ( status[i] == 0 ) {
                prime[primeSize++] = i;
            }
        }
    }

    void factorizeWithSieve(int n) {
        int sqrtn = sqrt(n);
        for ( int i = 0; i < primeSize && prime[i] <= sqrtn; i++ ) {
            if ( n % prime[i] == 0 ) {
                while ( n % prime[i] == 0 ) {
                    factors.push_back(prime[i]);
                    n /= prime[i];
                }
                sqrtn = sqrt(n);
            }
        }
        if ( n != 1 ) {
            factors.push_back(n);
        }
    }

    ll pollard_rho( ll n, ll c ) {
        ll y = 2, i = 1, k = 2, d;
        __int128 x = 2;
        while (true) {
            x = x * x % n + c;
            if (x >= n)	x -= n;
            d = __gcd((ll)x - y, n);
            if (d > 1) return d;
            if (++i == k) {
                y = x, k <<= 1;
            }
        }
        return n;
    }

    void factorize(ll n) {
        if (n == 1)
            return ;
        if (n < 1e+9) {
            factorizeWithSieve(n);
            return ;
        }
        if (millerRabin.isPrime(n)) {
            factors.push_back(n);
            return ;
        }
        ll d = n;
        for (int i = 2; d == n; i++) {
            d = pollard_rho(n, i);
        }
        factorize(d);
        factorize(n/d);
    }

    public:

    vector<ll> factors;

    PollardRho() {
        sieve();
    }

    void clear() {
        factors.clear();
    }

    vector<pair<ll,int>> getPrimeFactorization(ll n) {
        factorize(n);
        sort(factors.begin(), factors.end());

        vector<pair<ll,int>> res;
        for( int i = 0; i < factors.size(); i++ ) {
            ll p = factors[i];
            int cnt = 1;
            while ( i + 1 < factors.size() && factors[i+1] == p) {
                i++;
                cnt++;
            }
            res.push_back({p,cnt});
        }

        return res;
    }
}pollardRho;

/***************************/

int main() {
    ll n = 1e18;

    pollardRho.clear(); /// Don't forget to clear. Important for multi case.
    vector<pair<ll,int>> factors = pollardRho.getPrimeFactorization(n);
    for ( int i = 0; i < factors.size(); i++ ) {
        ll p = factors[i].first;
        ll a = factors[i].second;

        cout << p << ' ' << a << endl;

        /// p^a is factor of n
        /// Do your work here
    }
}
