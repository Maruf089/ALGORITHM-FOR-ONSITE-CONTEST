#include <bits/stdc++.h>
using namespace std;
#define ll long long int
int n,m,a,b,t,i,j,d,cs=0,counT=0,k,ans=0,l=0,sum1=0,sum=0,Max = -1,Min,num;
typedef  pair<int, int> iPair;
struct Graph
{
    int V, E;
    vector< pair<ll, iPair> > edges;

    // Constructor
    Graph(int V)
    {
        this->V = V;
    }

    // Utility function to add an edge
    void addEdge(int u, int v,ll w)
    {
        edges.push_back({w, {u, v}});
    }
    ll kruskalMST();
};

// To represent Disjoint Sets
struct DisjointSets
{
    int *parent, *rnk;
    int n;

    // Constructor.
    DisjointSets(int n)
    {
        // Allocate memory
        this->n = n;
        parent = new int[n+1];
        rnk = new int[n+1];

        // Initially, all vertices are in
        // different sets and have rank 0.
        for (int i = 0; i <= n; i++)
        {
            rnk[i] = 0;

            //every element is parent of itself
            parent[i] = i;
        }
    }

    // Find the parent of a node 'u'
    // Path Compression
    int find(int u)
    {
        /* Make the parent of the nodes in the path
           from u--> parent[u] point to parent[u] */
        if (u != parent[u])
            parent[u] = find(parent[u]);
        return parent[u];
    }

    // Union by rank
    void merge(int x, int y)
    {
        x = find(x), y = find(y);

        /* Make tree with smaller height
           a subtree of the other tree  */
        if (rnk[x] > rnk[y])
            parent[y] = x;
        else // If rnk[x] <= rnk[y]
            parent[x] = y;

        if (rnk[x] == rnk[y])
            rnk[y]++;
    }
};

ll Graph::kruskalMST()
{
    ll mst_wt = 0; // Initialize result
    // Sort edges in increasing order on basis of cost
    sort(edges.rbegin(), edges.rend());
    // Create disjoint sets
    DisjointSets ds(V);

    // Iterate through all sorted edges
    vector< pair<ll, iPair> >::iterator it;
    for (it=edges.begin(); it!=edges.end(); it++)
    {
        int u = it->second.first;
        int v = it->second.second;

        int set_u = ds.find(u);
        int set_v = ds.find(v);

        // Check if the selected edge is creating
        // a cycle or not (Cycle is created if u
        // and v belong to same set)
        if (set_u != set_v)
        {
            // Current edge will be in the MST
            // so print it
            //  cout << u << " - " << v << endl;

            // Update MST weight
            mst_wt += it->first;

            // Merge two sets
            ds.merge(set_u, set_v);
        }
    }

    return mst_wt;
}

int main()
{
    cin >> n >> d ;
    vector<vector<int>>point(n,vector<int>(d));
    Graph g(n);

    for(int i=0; i<n; i++)
    {
        for(int j=0; j<d; j++)
        {
            cin >> point[i][j];
        }
    }

    for(int mask=0; mask<(1<<d); mask++)
    {
        vector<pair<int,int>>e;
        for(int i=0; i<n; i++)
        {
            int sum = 0 ;
            for(int j=0; j<d; j++)
            {
                if((mask>>j)&1)
                    sum += point[i][j];
                else
                    sum -= point[i][j];
            }
            e.push_back({sum,i});
        }

        auto func = [&](int u,int v)
        {
            ll sum = 0;
            for(int i=0; i<d; i++)
                sum += abs(point[u][i]-point[v][i]);
            g.addEdge(u,v,sum);
        };


        int small = min_element(e.begin(),e.end())->second;
        int large = max_element(e.begin(),e.end())->second;
        for(int x=0 ; x<n; x++)
        {
            func(small,x);
            func(large,x);
        }
    }
    ll ans = g.kruskalMST();
    cout << ans << endl;

}
