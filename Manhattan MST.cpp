///*بِسْمِ اللَّهِ الرَّحْمَنِ الرَّحِيم*///

//#pragma GCC optimize("Ofast")
//#pragma GCC target("avx,avx2,fma")
//#pragma GCC optimization ("unroll-loops")

#include <bits/stdc++.h>
using namespace std;
#define ln '\n'
#define inp(x) scanf("%lld",&x)
#define inp2(a,b) scanf("%lld %lld",&a,&b)
#define sf scanf
#define pf printf

#define f0(i,b) for(int i=0;i<(b);i++)
#define f1(i,b) for(int i=1;i<=(b);i++)

#define MOD 1000000007
#define PI acos(-1)
#define ll long long int
typedef pair<int,int> pii;
typedef pair<ll, ll> pll;
#define pb push_back
#define F first
#define S second
#define sz size()
#define LCM(a,b) (a*(b/__gcd(a,b)))
//#define harmonic(n) 0.57721566490153286060651209+log(n)+1.0/(2*n)
#define SORT_UNIQUE(c) (sort(c.begin(),c.end()), c.resize(distance(c.begin(),unique(c.begin(),c.end()))))
#define all(v) v.begin(),v.end()
#define MEM(a, b) memset(a, b, sizeof(a))

#define fastio  ios_base::sync_with_stdio(false);
///Inline functions

inline bool EQ(double a, double b) { return fabs(a-b) < 1e-9; }
//inline bool isLeapYeall year) { return (year%400==0) | (year%4==0 && year%100!=0); }
inline void normal(ll &a) { a %= MOD; (a < 0) && (a += MOD); }
inline ll modMul(ll a, ll b) { a %= MOD, b %= MOD; normal(a), normal(b); return (a*b)%MOD; }
inline ll modAdd(ll a, ll b) { a %= MOD, b %= MOD; normal(a), normal(b); return (a+b)%MOD; }
inline ll modSub(ll a, ll b) { a %= MOD, b %= MOD; normal(a), normal(b); a -= b; normal(a); return a; }
inline ll modPow(ll b, ll p) { ll r = 1; while(p) { if(p&1) r = modMul(r, b); b = modMul(b, b); p >>= 1; } return r; }
inline ll modInverse(ll a) { return modPow(a, MOD-2); }
inline ll modDiv(ll a, ll b) { return modMul(a, modInverse(b)); }
inline bool isInside(pii p,ll n,ll m){ return (p.first>=0&&p.first<n&&p.second>=0&&p.second<m); }
inline bool isInside(pii p,ll n){ return (p.first>=0&&p.first<n&&p.second>=0&&p.second<n); }
inline bool isSquare(ll x){ ll s = sqrt(x);    return (s*s==x); }
inline bool isFib(ll x) { return isSquare(5*x*x+4)|| isSquare(5*x*x-4); }
inline bool isPowerOfTwo(ll x){ return ((1LL<<(ll)log2(x))==x); }


/// DEBUG --------------------------------------------------------------------------------->>>>>>

///**
template < typename F, typename S >
ostream& operator << ( ostream& os, const pair< F, S > & p )
{
    return os << "(" << p.first << ", " << p.second << ")";
}

template < typename T >
ostream &operator << ( ostream & os, const vector< T > &v )
{
    os << "{";
    for (auto it = v.begin(); it != v.end(); ++it)
    {
        if ( it != v.begin() )
            os << ", ";
        os << *it;
    }
    return os << "}";
}

template < typename T >
ostream &operator << ( ostream & os, const set< T > &v )
{
    os << "[";
    for (auto it = v.begin(); it != v.end(); ++it)
    {
        if ( it != v.begin())
            os << ", ";
        os << *it;
    }
    return os << "]";
}

template < typename F, typename S >
ostream &operator << ( ostream & os, const map< F, S > &v )
{
    os << "[";
    for (auto it = v.begin(); it != v.end(); ++it)
    {
        if ( it != v.begin() )
            os << ", ";
        os << it -> first << " = " << it -> second ;
    }
    return os << "]";
}

#define dbg(args...) do {cerr << #args << " : "; faltu(args); } while(0)

clock_t tStart = clock();
#define timeStamp dbg("Execution Time: ", (double)(clock() - tStart)/CLOCKS_PER_SEC)

void faltu ()
{
    cerr << endl;
}

template <typename T>
void faltu( T a[], int n )
{
    for (int i = 0; i < n; ++i)
        cerr << a[i] << ' ';
    cerr << endl;
}

template <typename T, typename ... hello>
void faltu( T arg, const hello &... rest)
{
    cerr << arg << ' ';
    faltu(rest...);
}

/// TEMPLATE ----------------------------------------------------------------------->>>>>>


struct func
{
//this is a sample overloading function for sorting stl
    bool operator()(pii const &a, pii const &b)
    {
        if(a.F==b.F)
            return (a.S<b.S);
        return (a.F<b.F);
    }
};
/*------------------------------Graph Moves----------------------------*/
//Rotation: S -> E -> N -> W
//const int fx[] = {0, +1, 0, -1};
//const int fy[] = {-1, 0, +1, 0};
//const int fx[]={+0,+0,+1,-1,-1,+1,-1,+1};   // Kings Move
//const int fy[]={-1,+1,+0,+0,+1,+1,-1,-1};  // Kings Move
//const int fx[]={-2, -2, -1, -1,  1,  1,  2,  2};  // Knights Move
//const int fy[]={-1,  1, -2,  2, -2,  2, -1,  1}; // Knights Move
/*---------------------------------------------------------------------*/

/// Bit Operations
/// count set bit C = (num * 0x200040008001ULL & 0x111111111111111ULL) % 0xf; /// 32 bit integer
/// inline bool checkBit(ll n, int i) { return n&(1LL<<i); }
/// inline ll setBit(ll n, int i) { return n|(1LL<<i); }
/// inline ll resetBit(ll n, int i) { return n&(~(1LL<<i)); }


/// ************************************** Code starts here ****************************************** */

/*
// Program showing a policy-based data structure.
#include <ext/pb_ds/assoc_container.hpp> // Common file
#include <ext/pb_ds/tree_policy.hpp>
#include <functional> // for less
using namespace __gnu_pbds;

// GNU link : https://goo.gl/WVDL6g
typedef tree<int, null_type, less_equal<int>, rb_tree_tag,
        tree_order_statistics_node_update>
        new_data_set;
// find_by_order(k) – ফাংশনটি kth ordered element এর একটা পয়েন্টার রিটার্ন করে। অর্থাৎ তুমি চাইলেই kth ইন্ডেক্সে কি আছে, সেটা জেনে ফেলতে পারছো!
// order_of_key(x) – ফাংশনটি x এলিমেন্টটা কোন পজিশনে আছে সেটা বলে দেয়।
*/


const int N = 200006;
int n;
vector<pii> g[N];
struct PT
{
  int x, y, id;
  bool operator < (const PT &p) const {
    return x == p.x ? y < p.y : x <  p.x;
  }
} p[N];

struct node
{
  int val, id;
} t[N];

struct DSU {
  int p[N];
  void init(int n) { for (int i = 1; i <= n; i++)  p[i] = i; }
  int find(int u) { return p[u] == u ? u : p[u] = find(p[u]); }
  void merge(int u, int v) { p[find(u)] = find(v); }
} dsu;

struct edge
{
  int u, v, w;
  bool operator < (const edge &p) const {return w < p.w;}
};
vector<edge> edges;

int Query(int x) {
  int r = 2e9 + 20;
  int id = -1;
  for (; x <= n; x += (x & -x))
    if (t[x].val < r) {
      r = t[x].val;
      id = t[x].id;
    }
  return id;
}

void modify(int x, int w, int id) {
  for (; x > 0; x -= (x & -x))
    if (t[x].val > w) {
      t[x].val = w;
      t[x].id = id;
    }
}

int dist(PT &a, PT &b) {
  return abs(a.x - b.x) + abs(a.y - b.y);
}

void add(int u, int v, int w) {
  edges.push_back({u, v, w});
}

ll Kruskal() {
  dsu.init(n);
  sort(edges.begin(), edges.end());
  ll ans = 0;
  for (edge e : edges) {
    int u = e.u, v = e.v, w = e.w;
    if (dsu.find(u) != dsu.find(v)) {
      ans += w;
      g[u].push_back({v, w});
      dsu.merge(u, v);
    }
  }
  return ans;
}

void Manhattan() {
  for (int i = 1; i <= n; i++)p[i].id = i;
  for (int dir = 1; dir <= 4; ++dir) {
    if (dir == 2 || dir == 4) {
      for (int i = 1; i <= n; i++)
        swap(p[i].x, p[i].y);
    }
    else if (dir == 3)
      for (int i = 1; i <= n; i++)
        p[i].x = -p[i].x;
    sort(p + 1, p + 1 + n);
    vector<int> v;
    static int a[N];
    for (int i = 1; i <= n; i++) {
      a[i] = p[i].y - p[i].x;
      v.push_back(a[i]);
    }
    sort(v.begin(), v.end());
    v.erase(unique(v.begin(), v.end()), v.end());
    for (int i = 1; i <= n; i++) {
      a[i] = lower_bound(v.begin(), v.end(), a[i]) - v.begin() + 1;
    }
    for (int i = 1; i <= n; i++) {
      t[i].val = 2e9 + 10;
      t[i].id = -1;
    }
    for (int i = n; i >= 1; i--) {
      int pos = Query(a[i]);
      if (pos != -1)
        add(p[i].id, p[pos].id, dist(p[i], p[pos]));
      modify(a[i], p[i].x + p[i].y, i);
    }
  }
}

void solve() {
  cin >> n ;
  for (int i = 1; i <= n; i++)
  {
      cin >> p[i].x >> p[i].y;
  }
  Manhattan();
  cout << Kruskal() << "\n";
  for (int u = 1; u <= n; u++) {
    for (auto &x : g[u])
      cout << u - 1 << " " << x.first - 1 << "\n";
  }
  return;
}

int main()
{
  fastio;
  solve();
}

