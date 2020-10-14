
#include<bits/stdc++.h>
using namespace std;

clock_t tStart = clock();
#define timeStamp cout << (double)(clock() - tStart)/CLOCKS_PER_SEC << endl;

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
