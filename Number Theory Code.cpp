

/// divisor count from 1 to N /// here MAX = 10^7
/// Extended Euclid , /// finds x and y , ax+by = gcd(a,b).
/// Bigmod
/// Find modinverse from 1 to N in O(n)
/// Sieve of Erathosnes
/// BITWISE SIEVE
/// CountDiv()
/// primeFactor()
/// Segmented Sieve ,  find number of primes inside the range of A and B inclusive. Here, A,B≤10^12 and  B–A≤10^6
/// Segmented Phi
/// Find All prime between 2 : N in O(n)
/// NOD , number of divisor of a number O(sqrt)
/// SOD , sum over the divisor function. O(sqrt)
/// Linear Diophantine Equation, Ax + By = C ( find x, y )
/// Simple Hyperbolic Diophantine Equation ,  Axy + Bx + Cy = D
/// Catalan number ( 1 1 2 5 14 42 132 429 1430 4862 )
/// seive phi..for Precalculate
/// single number phi in O(sqrt(n))
/// binomial coefficient( NcR in O(r) )
/// factorial inverse from 1 to MAXN in O(MAXN)
/// n'th fibonacci
/// prime count <= 10^11 /// Complexity for <= M , O(M^(3/4))
/// prime count /// Complexity sqrt(n)
/// count and print prime <= 5*10^8 ( 7104 ms )
/// Radix sort
/// Number base system
/// digits in factorial N!
/// prime factorization of N!
/// Leading digits of N!
/// Chinese Remainder Theorem , also another code Works for non-co-prime moduli




#include<bits/stdc++.h>
using namespace std;

#define NN 1000005
long total[1000005];
bool Isprime[NN];
int prime[NN];
int totalPrime;


void div()
{
    for(int i=1; i<MAX; i++) /// here MAX = 10^7
    {
        for(int j=i+i; j<MAX; j=j+i)
        {
           /// i is a divisor of j
           divisor[j]++;
        }
    }
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


/// x,y global
/// It finds two integers x and y such that, ax+by = gcd(a,b).
void extEuclid(int64 a, int64 b) {
	if (b == 0)
        { x = 1; y = 0; d = a; return; }
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

/// Find modinverse from 1 to N in O(n)
vector<int> inverseArray(int n, int m) {
    vector<int> modInverse(n + 1,0);
    modInverse[1] = 1;
    for(int i = 2; i <= n; i++) {
        modInverse[i] = (-(m/i) * modInverse[m % i]) % m + m;
    }
    return modInverse;
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
            {
                Set(prime[j>>5],j&31);
            }
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
int main()
{
	int i,j;
   PHI();
	int n;
	while(scanf("%d",&n),n)
	{
		printf("%d\n",phi[n]);
	    cout << (phi[n]/2)*n << endl; /// Sum of Co-prime Numbers of an Integer
	}
	return 0;
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

/// Segmented Phi
const int SIZE = (int)1e5+1;     //define size to be max(B-A+1)
ll r[SIZE];
ll phi[SIZE];
void segmented_phi(ll a, ll b){
    for(int i = 0;i < SIZE; i++){
        r[i] = a+i;
        phi[i] = a+i;
    }

    for(int i = 0;i < primes.size(); i++){
        int p = primes[i];
        ll j = p;
        if(j < a)    j = ((a+p-1)/p)*p;    //j = ceil(a/p)*p
        for(; j <= b; j+=p){    //j is the smallest multiple of p
            phi[j-a] /= p;
            phi[j-a] *= (p-1);
            while(r[j-a]%p == 0){
                r[j-a]/=p;
            }
        }
    }
    for(ll i = a; i<= b; i++){
        if(r[i-a] > 1){
            phi[i-a] /= r[i-a];
            phi[i-a] *= (r[i-a]-1);
        }
    }
    for(ll i = a;i <= b; i++){
      printf("%lld\n",phi[i-a]);
    }
    return;
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

/// number of divisor of a number O(sqrt)
int NOD ( int n ) {
    int sqrtn = sqrt ( n );
    int res = 1;
    for ( int i = 0; i < prime.size() && prime[i] <= sqrtn; i++ ) {
        if ( n % prime[i] == 0 ) {
            int p = 0; // Counter for power of prime
            while ( n % prime[i] == 0 ) {
                n /= prime[i];
                p++;
            }
            sqrtn = sqrt ( n );
            p++; // Increase it by one at end
            res *= p; // Multiply with answer
        }
    }
    if ( n != 1 ) {
        res *= 2; // Remaining prime has power p^1. So multiply with 2/
    }
    return res;
}


/// sum of divisor of a single number O(sqrt)
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
        res *= ( n + 1 ); /// Need to multiply (p^0+p^1)
    }
    return res;
}


/// The function CSOD(n) (cumulative SOD) of aninteger n, is defined as below:
///  CSOD(n)=  { from(i=1 to n) SOD(i) }
ll CSOD(int n)
{
    ll sum =0  ;
    for(i=2;i<=sqrt(n);i++)
    {
        j = n/i;
        sum += ((j*(j+1))-(i*(i-1)))/2;
        sum += (j-i)*i;
    }
    cout<< sum << endl;
}

/// the divisor summatory function is a function that is a sum over the divisor function. O(sqrt)
int SNOD( int n ) {
    int res = 0;
    int u = sqrt(n);
    for ( int i = 1; i <= u; i++ ) {
        res += ( n / i ) - i; //Step 1
    }
    res *= 2; //Step 2
    res += u; //Step 3
    return res;
}


/// Linear Diophantine Equation
/// Ax + By = C ( find x, y )
bool linearDiophantine ( int A, int B, int C, int *x, int *y ) {
    int g = gcd ( A, B );
    if ( C % g != 0 ) return false; //No Solution

    int a = A / g, b = B / g, c = C / g;
    ext_gcd( a, b, x, y ); //Solve ax + by = 1

    if ( g < 0 ) { //Make Sure gcd(a,b) = 1
        a *= -1; b *= -1; c *= -1;
    }

    *x *= c; *y *= c; //ax + by = c
    return true; //Solution Exists
}

int main () {
    int x, y, A = 2, B = 3, C = 5;
    bool res = linearDiophantine ( A, B, C, &x, &y );

    if ( res == false ) printf ( "No Solution\n" );
    else {
        printf ( "One Possible Solution (%d %d)\n", x, y );

        int g = gcd ( A, B );

        int k = 1; //Use different value of k to get different solutions
        printf ( "Another Possible Solution (%d %d)\n", x + k * ( B / g ), y - k * ( A / g ) );
    }

 return 0;
}



/// Simple Hyperbolic Diophantine Equation
/// Axy + Bx + Cy = D

bool isValidSolution ( int a, int b, int c, int p, int div ) {
    if ( ( ( div - c )% a ) != 0 ) return false; //x = (div - c) / a
    if ( ( (p-b*div) % (a*div) ) != 0 ) return false;// y = (p-b*div) /(a*div)
    return true;
}

int hyperbolicDiophantine ( int a, int b, int c, int d ) {
    int p = a * d + b * c;

    if ( p == 0 ) { //ad + bc = 0
        if ( -c % a == 0 ) return -1; //Infinite solutions (-c/a, k)
        else if ( -b % a == 0 ) return -1; //Infinite solutions (k, -b/a)
        else return 0; //No solution
    }
    else {
        int res = 0;

        //For each divisor of p
        int sqrtn = sqrt ( p ), div;
        for ( int i = 1; i <= sqrtn; i++ ) {
            if ( p % i == 0 ) { //i is a divisor
                //Check if divisors i,-i,p/i,-p/i produces valid solutions
                if ( isValidSolution(a,b,c,p,i) )res++;
                if ( isValidSolution(a,b,c,p,-i) )res++;
                if ( p / i != i ) { //Check whether p/i is different divisor than i
                    if ( isValidSolution(a,b,c,p,p/i) )res++;
                    if ( isValidSolution(a,b,c,p,-p/i) )res++;
                }
            }
        }

        return res; /// return number of solution exist
    }
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


/// seive phi..
ll phi[MAX];
void seivePHI()
{
    /// phi[1] = 1;
    for(i = 2; i < MAX; i++)
    {
        if(phi[i] == 0)
        {
            phi[i] = i - 1;
            for(j = i*2; j < MAX; j += i)
            {
                if(phi[j] == 0) phi[j] = j;
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


long long NCR(int n, int r) { /// O(r)
    if(r > n - r) r = n - r; // because C(n, r) == C(n, n - r)
    long long ans = 1;
    for(i = 1; i <= r; i++) {
        ans *= n - r + i; /// or ans *= n-i+1
        ans /= i;
    }
    return ans;
}
/// if N and R large but MOD value small we can find NcR in  O(MOD)
ll fact_inv() /// O(n)
{
    fact[0 = inv[0] = 1;
    for(i=1;i<maxn;i++) fact[i] = modMul(fact[i-1],i);
    inv[maxn-1] = modPow(fact[maxn-1],MOD-2);
    for(i=maxn-2;i>=1;i--)inv[i] = modMul(inv[i+1],i+1);
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



/// count and print prime <= 5*10^8 ( 7104 ms )
bitset<500000009>sieve;
vector<int>result;
int main()
{
    long long N = 500000000, A, B;
    int number_of_prime = -1;
    for (int i = 2; i <= N; ++i)
    {
        if (sieve[i] == 0)
        {
            ++number_of_prime;
            for (int j = i + i; j <= N; j+=i)
                sieve[j] = 1;
            result.push_back(i);
        }
    }
    cout << result.size();
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



/// Number base system
// A list of symbol. Depending on base and number system, this list can be different.
char symbol[] = {'0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F'};
string decimalToBase ( int x, int base ) {
    string res = "";

    while ( x ) {
        int r = x % base; // Find the last digit
        res = res + symbol[r]; // Change the integer value to symbol and append to res
        x /= base; // Remove the last digit
    }
    if ( res == "" ) res = symbol[0]; // If res is empty, that means x is 0.
    reverse ( res.begin(), res.end()); // We found the digits in reverse order.
    return res;
}


/// digits in factorial N
int factorialDigitExtended ( int n, int base ) {
    double x = 0;
    for ( int i = 1; i <= n; i++ ) {
        x += log10 ( i ) / log10(base); // Base Conversion
    }
    int res = x + 1 + eps;
    return res;
}

/// prime factorization of N!
void factFactorize ( int n ) {
    for ( int i = 0; i < prime.size() && prime[i] <= n; i++ ) {
        int x = n;
        int freq = 0;

        while ( x / prime[i] ) {
            freq += x / prime[i];
            x = x / prime[i];
        }

        printf ( "%d^%d\n", prime[i], freq );
    }
}

/// Leading digits of N!

const double eps = 1e-9;
/// Find the first K digits of N!
int leadingDigitFact ( int n, int k ) {
    double fact = 0;

    /// Find log(N!)
    for ( int i = 1; i <= n; i++ ) {
        fact += log10 ( i );
    }

    /// Find the value of q
    double q = fact - floor ( fact+eps );

    double B = pow ( 10, q );

    /// Shift decimal point k-1 \times
    for ( int i = 0; i < k - 1; i++ ) {
        B *= 10;
    }

    /// Don't forget to floor it
    return floor(B+eps);
}


/// Chinese Remainder Theorem
/** Return {-1,-1} if invalid input.
    Otherwise, returns {x,L}, where x is the solution unique to mod L
*/
pair<int, int> chinese_remainder_theorem( vector<int> A, vector<int> M ) {
    if(A.size() != M.size()) return {-1,-1}; /** Invalid input*/

    int n = A.size();

    int a1 = A[0];
    int m1 = M[0];
    /** Initially x = a_0 (mod m_0)*/

    /** Merge the solution with remaining equations */
    for ( int i = 1; i < n; i++ ) {
        int a2 = A[i];
        int m2 = M[i];

        /** Merge the two equations*/
        int p, q;
        ext_gcd(m1, m2, &p, &q);

        /** We need to be careful about overflow, but I did not bother about overflow here to keep the code simple.*/
        int x = (a1*m2*q + a2*m1*p) % (m1*m2);

        /** Merged equation*/
        a1 = x;
        m1 = m1 * m2;
    }
    if (a1 < 0) a1 += m1; /** Result is not suppose to be negative*/
    return {a1, m1};
}


/** Works for non-coprime moduli.
 Returns {-1,-1} if solution does not exist or input is invalid.
 Otherwise, returns {x,L}, where x is the solution unique to mod L
*/
pair<int, int> chinese_remainder_theorem( vector<int> A, vector<int> M ) {
    if(A.size() != M.size()) return {-1,-1}; /** Invalid input*/

    int n = A.size();

    int a1 = A[0];
    int m1 = M[0];
    /** Initially x = a_0 (mod m_0)*/

    /** Merge the solution with remaining equations */
    for ( int i = 1; i < n; i++ ) {
        int a2 = A[i];
        int m2 = M[i];

        int g = __gcd(m1, m2);
        if ( a1 % g != a2 % g ) return {-1,-1}; /** No solution exists*/

        /** Merge the two equations*/
        int p, q;
        ext_gcd(m1/g, m2/g, &p, &q);

        int mod = m1 / g * m2; /** LCM of m1 and m2*/

        /** We need to be careful about overflow, but I did not bother about overflow here to keep the code simple.*/
        int x = (a1*(m2/g)*q + a2*(m1/g)*p) % mod;

        /** Merged equation*/
        a1 = x;
        if (a1 < 0) a1 += mod; /** Result is not suppose to be negative*/
        m1 = mod;
    }
    return {a1, m1};
}

