#include<bits/stdc++.h>
using namespace std;
inline long bit_block_count(long x)
{
   return (x & 1) + __builtin_popcountll((x^(x>>1)) ) / 2;
}
inline long bit_block_greater_than_2_count(long x)
{
    return bit_block_count( x & ( (x<<1) & (x>>1) ) );
}

int __builtin_ffs (unsigned int x)
/// Returns one plus the index of the least significant 1-bit of x,
/// or if x is zero, returns zero.
int __builtin_clz (unsigned int x)
/// Returns the number of leading 0-bits in x, starting at the
/// most significant bit position. If x is 0, the result is undefined.
int __builtin_ctz (unsigned int x)
/// Returns the number of trailing 0-bits in x, starting at the
/// least significant bit position. If x is 0, the result is undefined.
int __builtin_popcount (unsigned int x)
/// Returns the number of 1-bits in x.
int __builtin_parity (unsigned int x)
/// Returns the parity of x, i.e. the number of 1-bits in x modulo 2.

struct custom_hash {
    static uint64_t splitmix64(uint64_t x) {
        // http://xorshift.di.unimi.it/splitmix64.c
        x += 0x9e3779b97f4a7c15;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9;
        x = (x ^ (x >> 27)) * 0x94d049bb133111eb;
        return x ^ (x >> 31);
    }

    size_t operator()(uint64_t x) const {
        static const uint64_t FIXED_RANDOM = chrono::steady_clock::now().time_since_epoch().count();
        return splitmix64(x + FIXED_RANDOM);
    }
};

unsigned int nextPowerOf2(unsigned int n)
{
    n--;
    n |= n >> 1;
    n |= n >> 2;
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16;
    n++;
    return n;
}

int highestPowerof2(int n)
{
   int p = (int)log2(n);
   return (int)pow(2, p);
}

int msb(int x) {
    union { double a; int b[2]; };
    a = x;
    return (b[1] >> 20) - 1023;
}

int highestOneBit(int i) {
    // HD, Figure 3-1
    i |= (i >>  1);
    i |= (i >>  2);
    i |= (i >>  4);
    i |= (i >>  8);
    i |= (i >> 16);
    return i - ((unsigned)i >> 1);
}
static inline ulong bit_swap(ulong a, ulong k1, ulong k2)
 /// Return a with bits at positions [k1] and [k2] swapped.
 /// k1==k2 is allowed (a is unchanged then)
 {
 ulong x = ((a>>k1) ^ (a>>k2)) & 1; // one if bits differ
 a ^= (x<<k2); // change if bits differ
 a ^= (x<<k1); // change if bits differ
 return a;
}

int getFirstSetBitPos(int n){ return ffs(n);}
int main()
{
cout << bit_block_count(37) << " " << bit_block_greater_equal_2_count(79);
}


/// Bitwise Sieve
#define M 100000100
int flag[M/32];
int cnt;
int prime[5761482];
unsigned ans;
unsigned store[5761482];

void prime_gen()
{
    int add,x=0;
    prime[x++]=2;
    for(int i = 4; i<M; i+=2)
        flag[i/32]=_set(flag[i/32],i%32);
    int sq = sqrt(M);
    for(int i = 3; i<M; i+=2)
    {
        if(check(flag[i/32],i%32)==0)
        {
            prime[x++]=i;
            if(sq>=i)
            {
                add = i*2;
                for(int j = i*i; j<M; j+=add)
                    flag[j/32]=_set(flag[j/32],j%32);
            }
        }

    }
    cnt=x;
}


/*
#include <ext/pb_ds/assoc_container.hpp> // Common file
#include <ext/pb_ds/tree_policy.hpp> // Including */

//using namespace __gnu_pbds;
/*
template<typename T>
using ordered_set = tree<T, null_type, less<T>, rb_tree_tag, tree_order_statistics_node_update>;
template<typename F, typename S>
using ordered_map = tree<F, S, less<F>, rb_tree_tag, tree_order_statistics_node_update>;

// find_by_order(k) – ফাংশনটি kth ordered element এর একটা পয়েন্টার রিটার্ন করে। অর্থাৎ তুমি চাইলেই kth ইন্ডেক্সে কি আছে, সেটা জেনে ফেলতে পারছো!
// order_of_key(x) – ফাংশনটি x এলিমਠƨ্টটޠকোন পজিশনে আছে সেটা বলে দেয়।

*/
