#include <bits/stdc++.h>
using namespace std;

void counting_sort(vector<int>&Order,vector<int>&Class)
{
    int n = Order.size();
    vector<int>cnt(n);
    for(auto x : Class)
        cnt[x]++;

    vector<int>Order_new(n);
    vector<int>pos(n);
    pos[0] = 0;
    for(int i=1; i<n; i++)
        pos[i] = pos[i-1] + cnt[i-1];

    for(auto it : Order)
    {
        int i = Class[it];
        Order_new[pos[i]] = it;
        pos[i]++;
    }
    Order = Order_new;
}

int main()
{
    string s;
    cin >> s ;
    s += "$";
    int n = s.size();
    vector<int>Order(n), Class(n);

    {
        /// k = 0
        vector<pair<char,int>>A(n);
        for(int i=0; i<n; i++)
            A[i] = {s[i],i};

        sort(A.begin(),A.end());
        for(int i=0; i<n; i++)
            Order[i] = A[i].second;

        Class[Order[0]] = 0 ;
        for(int i=1; i<n; i++)
        {
            if(A[i].first==A[i-1].first)
                Class[Order[i]] = Class[Order[i-1]];
            else
                Class[Order[i]] = Class[Order[i-1]] + 1;
        }
    }

    {
        int k = 0;
        while( (1<<k) < n )
        {
            /// k -> k+1
            for(int i=0; i<n; i++)
                Order[i] = (Order[i] -  (1<<k) + n) % n;

            counting_sort(Order,Class);

            vector<int>new_Class(n);
            new_Class[Order[0]] = 0 ;
            for(int i=1; i<n; i++)
            {
                pair<int,int>prev = {Class[Order[i-1]],Class[(Order[i-1]+(1<<k))%n]};
                pair<int,int>now = {Class[Order[i]],Class[(Order[i]+(1<<k))%n]};
                if(prev==now)
                    new_Class[Order[i]] = new_Class[Order[i-1]];
                else
                    new_Class[Order[i]] = new_Class[Order[i-1]] + 1;
            }
            Class = new_Class;
            k++;
        }
    }

    vector<int>lcp(n);
    int k = 0;
    for(int i=0;i<n-1;i++)
    {
        int pi = Class[i];
        int j = Order[pi-1];
        /// lcp[i] = lcp(s[i...] , s[j....])
        while(s[i+k]==s[j+k])
            k++;
        lcp[pi] = k;
        k = max(0,k-1);
    }


    for(int i=0; i<n; i++)
    {
        cout << lcp[i] << ' ' << Order[i] <<  ' ' << s.substr(Order[i],n-Order[i]) << '\n';
    }

}
