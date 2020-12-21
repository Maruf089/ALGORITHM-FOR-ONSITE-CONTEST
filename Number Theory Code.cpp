#include<bits/stdc++.h>
using namespace std;


ll fact[N];
void factorial () {
    fact[0] = 1LL;
    for (ll i = 1; i < N; ++i) fact[i] = modMul(fact[i - 1], i);
}

ll NcR (ll n, ll r) {
    if (r > n) return 0LL;
    return modDiv(fact[n], modMul(fact[n-r], fact[r]));
}


#define NN 1000005
long total[1000005];
bool Isprime[NN];
int prime[NN];
int totalPrime;


void extEuclid(int64 a, int64 b) {
	if (b == 0) { x = 1; y = 0; d = a; return; }
	extEuclid(b, a % b);
	x = x - (a / b) * y;
	swap(x, y);
}

int bigmod(int n,int p,int m)
{
    if(p==0)
        return 1;
    int k = bigmod(n,p/2,m);

    k = (k*k)% m;
    if(p%2==1)
        k= (n*k) % m ;

    return k;
}

void sv()
{
    int t = sqrt(NN);
    Isprime[0]=Isprime[1]=1;

    for(int i = 4; i <= NN; i+= 2)
    Isprime[i] = true;

    for( int i=3; i<=t; i += 2 )
    {
        if( !Isprime[i] )
        {
            for( int j=i*i; j<NN; j+= i+i)
                Isprime[j] = true;
        }
    }
    totalPrime = 0;
    prime[totalPrime++] = 2;

    for(int i=3; i<NN; i+=2)
        if( !Isprime[i] )
            prime[totalPrime++] = i;
}



/// BITWISE SIEVE

bool Check(int N,int pos){ return (bool)(N & (1<<pos));}
void Set(int &N,int pos){ N=N | (1<<pos);}
const int MAX = 1e8+5;
int prime[(MAX>>5)+2];
vector < int > primes;
void sieve()
{
    int lim = sqrt(MAX);
    for(int i=3; i<=lim; i+=2)
    {
        if(!Check(prime[i>>5],i&31))
        {
            for(int j=i*i; j<=MAX; j+=(i<<1))
                Set(prime[j>>5],j&31);
        }
    }
    primes.push_back(2);
    for(int i=3; i<=MAX; i+=2)
    {
        if(!Check(prime[i>>5],i&31))
            primes.push_back(i);
    }
}


void CountDiv()
{
    for( int i=2;i<NN;i++ )
    {
        int N = i;
        int t = sqrt(N);
        int res = 1;
        for( int j=0;;j++ )
        {
            if( t<prime[j] ) break;
            int cnt = 1;
            while( N%prime[j]==0 )
            {
                N/=prime[j] ;
                cnt++;
            }
            t = sqrt(N);
            res*=cnt;
        }
        if( N>1 ) res*=2;
        d[i] = res;
    }
}


#define N (int)1e7+9

vector<int>prime ;
bool vis[N] ;
long phi[N] ;

void PHI(){
    phi[1] = 1 ;
    for(int i = 2 ; i < N ; ++i){
        if(!vis[i]){
            prime.push_back(i) ;
            phi[i] = i - 1 ;
        }
        for(int j = 0 ; j < prime.size() and i * prime[j] < N ; ++j){
            vis[i * prime[j]] = 1 ;
            if(i % prime[j] == 0){
                phi[i * prime[j]] = phi[i] * prime[j] ;
                break ;
            }
            else {
                phi[i * prime[j]] = phi[i] * phi[prime[j]] ;
            }
        }
    }

}




void primeFactor()
{
    for( int i=2;i<NN;i++ )
    {
        int N = i;
        int t = sqrt(N);
        for( int j=0;;j++ )
        {

            if( t < prime[j] ) break;
            bool hasFactor = false;

            while( N % prime[j] == 0 )
            {
                N /= prime[j] ;
                hasFactor = true;
            }
            t = sqrt(N);
            if( hasFactor )
                factor[i].push_back(prime[j]);
        }
        if( N>1 )
            factor[i].push_back(N);
    }
}


/// Given two integers A and B, find number of primes inside the range of A and B inclusive. Here, 1≤A≤B≤10^12 and  B–A≤10^5

int arr[SIZE];
int segmentedSieve ( int a, int b ) {
    if ( a == 1 )
        a++;
    int sqrtn = sqrt ( b );
    memset ( arr, 0, sizeof arr ); // Make all index of arr 0.
    for ( int i = 0; i < prime.size() && prime[i] <= sqrtn; i++ ) {
        int p = prime[i];
        int j = p * p;
        // If j is smaller than a, then shift it inside of segment [a,b]
        if ( j < a ) j = ( ( a + p - 1 ) / p ) * p;

        for ( ; j <= b; j += p )
            arr[j-a] = 1; // mark them as not prime
    }
    int res = 0;
    for ( int i = a; i <= b; i++ ) {
        // If it is not marked, then it is a prime
        if ( arr[i-a] == 0 )
            res++;
    }
    return res;
}


/// Find All prime between 2 : N in O(n)
/*

const int N = 10000000;
int lp[N+1];
vector<int> pr;

for (int i=2; i<=N; ++i) {
    if (lp[i] == 0) {
        lp[i] = i;
        pr.push_back (i);
    }
    for (int j=0; j<(int)pr.size() && pr[j]<=lp[i] && i*pr[j]<=N; ++j)
        lp[i * pr[j]] = pr[j];
}
*/


int SOD( int n ) {
    int res = 1;
    int sqrtn = sqrt ( n );
    for ( int i = 0; i < prime.size() && prime[i] <= sqrtn; i++ ) {
        if ( n % prime[i] == 0 ) {
            int tempSum = 1; // Contains value of (p^0+p^1+...p^a)
            int p = 1;
            while ( n % prime[i] == 0 ) {
                n /= prime[i];
                p *= prime[i];
                tempSum += p;
            }
            sqrtn = sqrt ( n );
            res *= tempSum;
        }
    }
    if ( n != 1 ) {
        res *= ( n + 1 ); // Need to multiply (p^0+p^1)
    }
    return res;
}

int catalan[MAX];
void init() {
    catalan[0] = catalan[1] = 1;
    for (int i=2; i<=n; i++) {
        catalan[i] = 0;
        for (int j=0; j < i; j++) {
            catalan[i] += (catalan[j] * catalan[i-j-1]) % MOD;
            if (catalan[i] >= MOD) {
                catalan[i] -= MOD;
            }
        }
    }
}

int main()
{
	int i,j;
   PHI();
	int n;
	while(scanf("%d",&n),n)
	{
		printf("%d\n",phi[n]);
	//	cout << (phi[n]/2)*n << endl;
	}
	return 0;
}
/// seive phi..
ll phi[MAX];
void seivePHI()
{
    for(i = 2; i < MAX; i++)
    {
        if(phi[i] == 0)
        {
            phi[i] = i - 1;
            for(j = i*2; j < MAX; j += i)
            {
                if(phi[j] == 0)
                    phi[j] = j;
                phi[j] /= i;
                phi[j] *= (i-1);
            }
        }
    }
}

ll po(ll x, ll y)
{
    ans = 1;
    while(y--)
        ans *= x;
    return ans;
}
ll prime(ll a)
{
    for(ll i = 1; i*i <= a; i++)
    {
        if(a%i == 0)
            return 1;
    }
}
ll phi(ll n)
{
    ll i,mul = 1, holder, fre = 0;
    if(prime(n) == 0)
        mul = n - 1;
    else
    {
        for(i = 2; i*i <= n; i++)
        {
            if(n%i == 0)
            {
                while(n%i == 0)
                {
                    n = n/i;
                    holder = i;
                    fre++;
                }
                mul *= (po(holder, fre-1)*(holder - 1));
                fre = 0;
            }
        }
        if(n != 1)
        {
            mul *= (n-1);
        }
    }
    return mul;
}


int bc(int n, int k){
    int i, j;
    for (i = 0; i <= n; i++){
        for (j = 0; j <=  k; j++){
            if (j == 0 || j == i)
                DP[i][j] = 1;

            else
                DP[i][j] = (DP[i-1][j-1]%M + DP[i-1][j]%M)%M;
        }
    }
}


/// n'th fibonacci

ll fib(ll n)
{
    double phi = (1 + sqrt(5)) / 2;
    return round(pow(phi, n) / sqrt(5));
}


/// prime count <= 10^11

pi(n) = phi(n, sqrt(n)) + pi(sqrt(n)) - 1

where
   pi(n) = number of prime below N
   phi(n, m) = number of number below N that is not divisible by any prime below m.

#define rep(i,l,r)for(ll i=(l);i<(r);i++)
#define PLIMIT 1000000
int prime[PLIMIT+10];
void makep(int n){prime[0]=prime[1]=1;for(int i=2;i*i<=n;i++)if(!prime[i])for(int j=i*i;j<=n;j+=i)prime[j]=1;}

/// Complexity for <= M , O(M^(3/4))
ll*primecnt(ll M){
	ll SqM=sqrt(M);
	ll Len=SqM+M/(SqM+1)+1;
	ll*d=(ll*)malloc((Len+10)*sizeof(ll));
	makep(SqM+1);
	rep(i,1,SqM+1)d[i]=(M/i)-1;
	rep(i,SqM+1,Len)d[i]=(Len-i)-1;
	rep(i,0,SqM+1)if(!prime[i]){
		rep(n,1,Len){
			ll t=(n<=SqM?M/n:(Len-n))/i;
			if(t<i)break;
			d[n]-=(d[t>SqM?M/t:Len-t]-d[Len-(i-1)]);
		}
	}
	return d;
}

int main(){
	ll n;
	scanf("%lld",&n);
	printf("%lld\n",primecnt(n)[1]);
}


/// Complexity sqrt(n)

using i64 = long long;
int isqrt(i64 n) {
  return sqrtl(n);
}
__attribute__((target("avx"), optimize("O3", "unroll-loops")))
i64 prime_pi(const i64 N) {
  if (N <= 1) return 0;
  const int v = isqrt(N);
  vector<int> smalls(v + 1); for (int i = 2; i <= v; ++i) smalls[i] = (i + 1) / 2;
  int s = (v + 1) / 2;
  vector<int> roughs(s); for (int i = 0; i < s; ++i) roughs[i] = 2 * i + 1;
  vector<i64> larges(s); for (int i = 0; i < s; ++i) larges[i] = (N / (2 * i + 1) + 1) / 2;
  vector<bool> skip(v + 1);
  auto div = [] (i64 N, i64 d) -> int { return double(N) / d; };

  int pc = 0;
  for (int p = 3; p <= v; ++p) if (smalls[p] > smalls[p - 1]) {
    int q = p * p;
    pc++;
    if (i64(q) * q > N) break;
    skip[p] = true;
    for (int i = q; i <= v; i += 2 * p) skip[i] = true;
    int ns = 0;
    for (int k = 0; k < s; ++k) {
      int i = roughs[k];
      if (skip[i]) continue;
      i64 d = i64(i) * p;
      larges[ns] = larges[k] - (d <= v ? larges[smalls[d] - pc] : smalls[div(N, d)]) + pc;
      roughs[ns++] = i;
    }
    s = ns;
    for (int j = v / p; j >= p; --j) {
      int c = smalls[j] - pc;
      for (int i = j * p, e = min(i + p, v + 1); i < e; ++i) smalls[i] -= c;
    }
  }

  for (int k = 1; k < s; ++k) {
    const i64 M = N / roughs[k];
    i64 s = larges[k] - (pc + k - 1);
    for (int l = 1; l < k; ++l) {
      int p = roughs[l];
      if (i64(p) * p > M) break;
      s -= smalls[div(M, p)] - (pc + l - 1);
    }
    larges[0] -= s;
  }
  return larges[0];
}

int main(){
  i64 N;
  while (~scanf("%lld", &N)) {
    printf("%lld\n", prime_pi(N));
  }
}


/// Radix sort

int getMax(int arr[], int n)
{
    int mx = arr[0];
    for (int i = 1; i < n; i++)
        if (arr[i] > mx)
            mx = arr[i];
    return mx;
}

// A function to do counting sort of arr[] according to
// the digit represented by exp.
void countSort(int arr[], int n, int exp)
{
    int output[n]; // output array
    int i, count[10] = { 0 };

    // Store count of occurrences in count[]
    for (i = 0; i < n; i++)
        count[(arr[i] / exp) % 10]++;

    // Change count[i] so that count[i] now contains actual
    //  position of this digit in output[]
    for (i = 1; i < 10; i++)
        count[i] += count[i - 1];

    // Build the output array
    for (i = n - 1; i >= 0; i--) {
        output[count[(arr[i] / exp) % 10] - 1] = arr[i];
        count[(arr[i] / exp) % 10]--;
    }

    // Copy the output array to arr[], so that arr[] now
    // contains sorted numbers according to current digit
    for (i = 0; i < n; i++)
        arr[i] = output[i];
}

// The main function to that sorts arr[] of size n using
// Radix Sort
void radixsort(int arr[], int n)
{
    // Find the maximum number to know number of digits
    int m = getMax(arr, n);

    // Do counting sort for every digit. Note that instead
    // of passing digit number, exp is passed. exp is 10^i
    // where i is current digit number
    for (int exp = 1; m / exp > 0; exp *= 10)
        countSort(arr, n, exp);
}
