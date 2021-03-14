#define MAX_VERTEX 1000  // Maximum number of nodes a graph can have.
vector<int>G[MAX_VERTEX];// This basically is an array of Vector,
int color[MAX_VERTEX];
void dfs(int current_vertex)
{
    if(color[current_vertex])  // If the current vertex we are on is visited, just return to previously visited vertex.
    {
        return;
    }
    color[current_vertex]=1; // else color the current vertex.
    for(int i=0; i<G[current_vertex].size(); i++) // After coloring, just check it's neighbour,
    {
        // if they are not visited,visit them.
        if(!color[G[current_vertex][i]])
        {
            dfs(G[current_vertex][i]);
        }
    }
}


/// Connected Components
vector <int> adj[10];
bool visited[10];

void dfs(int s)
{
    visited[s] = true;
    for(int i = 0; i < adj[s].size(); ++i)
    {
        if(visited[adj[s][i]] == false)
            dfs(adj[s][i]);
    }
}
void initialize()
{
    for(int i = 0; i < 10; ++i)
        visited[i] = false;
}
int main()
{
    int nodes, edges, x, y, connectedComponents = 0;
    cin >> nodes;                       //Number of nodes
    cin >> edges;                       //Number of edges
    for(int i = 0; i < edges; ++i)
    {
        cin >> x >> y;
        //Undirected Graph
        adj[x].push_back(y);                   //Edge from vertex x to vertex y
        adj[y].push_back(x);                   //Edge from vertex y to vertex x
    }
    initialize();                           //Initialize all nodes as not visited
    for(int i = 1; i <= nodes; ++i)
    {
        if(visited[i] == false)
        {
            dfs(i);
            connectedComponents++;
        }
    }
    cout << "Number of connected components: " << connectedComponents << endl;
    return 0;
}






const int sz=10001;
vector<pair<int,int> > a[sz];
int dis[sz];
bool vis[sz]= {0};
int parent[sz];

void printPath(int j)
{
    if (parent[j] == - 1)
        return;
    printPath(parent[j]);
    printf("%d ", j);
}

void Dijkstra(int source, int n)
{
    for(int i=0; i<sz; i++)
        dis[i]=INF;
    parent[0] = -1 ;

    ///Custom Comparator for Determining priority for priority queue (shortest edge comes first)
    class prioritize
    {
    public:
        bool operator ()(pair<int, int>&p1,pair<int, int>&p2)
        {
            return p1.second>p2.second;
        }
    };
    priority_queue<pair<int,int>,vector<pair<int,int> >, prioritize> pq;  //Priority queue to store vertex,weight pairs


    pq.push(make_pair(source,dis[source]=0));
    while(!pq.empty())
    {
        pair<int, int> curr= pq.top(); //Current vertex. The shortest distance for this has been found
        pq.pop();
        int cv=curr.first,cw=curr.second; ///'cw' the final shortest distance for this vertex
        if(vis[cv]) ///If the vertex is already visited, no point in exploring adjacent vertices
            continue;
        vis[cv]=true;
        for(int i=0; i<a[cv].size(); i++)
        {
            int v = a[cv][i].first ;
            if(!vis[a[cv][i].first] && a[cv][i].second+cw<dis[a[cv][i].first]) //If this node is not visited and the current parent node distance+distance from there to this node is shorted than the initial distace set to this node, update it
            {
                pq.push(make_pair(a[cv][i].first,(dis[a[cv][i].first]=a[cv][i].second+cw))); //Set the new distance and add to priority queue
                parent[v] = a[cv][i].second;
            }
        }
    }
}

int main() ///Driver Function for Dijkstra SSSP
{
    int n,m,x,y,w;//Number of vertices and edges
    //cout<<"Enter number of vertices and edges in the graph\n";
    cin>>n>>m;
    for(int i=0; i<m; i++) //Building Graph
    {
        cin>>x>>y>>w; //Vertex1, Vertex2, weight of edge
        a[x].push_back(make_pair(y,w));
      //  a[y].push_back(make_pair(x,w));
    }
    //cout<<"Enter source for Dijkstra's SSSP algorithm\n";
    int source;
    cin>>source;
    Dijkstra(source,n);//SSSP from source (Also passing number of vertices as parameter)
    cout<<"Source is: "<<source<<". The shortest distance to every other vertex from here is: \n";
    for(int i=0; i<=n; i++) //Printing final shortest distances from source
    {
        cout<<"Vertex: "<<i<<" , Distance: ";
        dis[i]!=INF? cout<<dis[i]<<"\n" : cout<<"-1\n";
    }
}





void PrintNegativeCycle(vector< pair<int, int> > shortestDistances, int vertex, int startVertex)
{
    if (vertex == startVertex)
    {
        printf("%d ", vertex);
    }
    else if (shortestDistances[vertex].second == 0)
    {
        PrintNegativeCycle(shortestDistances, startVertex, startVertex);
        printf("%d ", vertex);
    }
    else
    {
        PrintNegativeCycle(shortestDistances, shortestDistances[vertex].second, startVertex);
        printf("%d ", vertex);
    }
}

// Bellman-Ford Algorithm which takes the Adjacency List, starting vertex,
// and an empty shortestDistances vector as input. It applies the algorithm
// and keeps filling values into shortestDistances which is a reference
// parameter. It returns true if there are no negative edges, and vice-versa.
int bellmanFord( vector< list< pair<int, int> > > adjacencyList,int vertices,int startVertex,vector< pair<int, int> > & shortestDistances)
{
    list< pair<int, int> >::iterator traverse;
    int i, j, k;

    // Initialisation
    for (i = 0; i <= vertices; ++i)
    {
        shortestDistances[i].first = INT_MAX;
        shortestDistances[i].second = -1;
    }

    // Setting distance to source = 0
    shortestDistances[startVertex].first = 0;
    shortestDistances[startVertex].second = 0;

    // The Algorithm that computes Shortest Distances
    for (i = 1; i <= vertices - 1; ++i)      // Runs 'vertices - 1' times = O(|V|)
    {
        for (j = 1; j <= vertices; ++j)      // Runs as many times as the edges = O(|E|)
        {
            // The code ahead basically explores the whole of Adjcency List,
            // covering one edge once, so the runtime of the code in this
            // block is O(|E|)

            traverse = adjacencyList[j].begin();

            while (traverse != adjacencyList[j].end())
            {
                if (shortestDistances[j].first == INT_MAX)
                {
                    // Important...!
                    //traverse = traverse->next;
                    ++traverse;
                    continue;
                }

                // Checking for Relaxation
                if ((*traverse).second + shortestDistances[j].first <
                        shortestDistances[(*traverse).first].first)
                {
                    // Relaxation
                    shortestDistances[(*traverse).first].first = (*traverse).second
                            + shortestDistances[j].first;
                    shortestDistances[(*traverse).first].second = j;
                }

                ++traverse;
            }
        }
    }

    // Checking for Negative Cycles
    for (j = 1; j <= vertices; ++j)
    {
        traverse = adjacencyList[j].begin();

        while (traverse != adjacencyList[j].end())
        {
            // Checking for further relaxation
            if ((*traverse).second + shortestDistances[j].first <
                    shortestDistances[(*traverse).first].first)
            {
                // Negative Cycle found as further relaxation is possible
                return j;
            }

            ++traverse;
        }
    }

    return -1;
}

int main()
{
    int vertices, edges, v1, v2, weight;

    printf("Enter the Number of Vertices -\n");
    scanf("%d", &vertices);

    printf("Enter the Number of Edges -\n");
    scanf("%d", &edges);

    // Adjacency List is a vector of list.
    // Where each element is a pair<int, int>
    // pair.first -> the edge's destination
    // pair.second -> edge's weight
    vector< list< pair<int, int> > > adjacencyList(vertices + 1);

    printf("Enter the Edges V1 -> V2, of weight W\n");

    for (int i = 1; i <= edges; ++i)
    {
        scanf("%d%d%d", &v1, &v2, &weight);

        // Adding Edge to the Directed Graph
        adjacencyList[v1].push_back(make_pair(v2, weight));
    }

    printf("\nThe Adjacency List-\n");
    // Printing Adjacency List
    for (int i = 1; i < adjacencyList.size(); ++i)
    {
        printf("adjacencyList[%d] ", i);

        list< pair<int, int> >::iterator itr = adjacencyList[i].begin();

        while (itr != adjacencyList[i].end())
        {
            printf(" -> %d(%d)", (*itr).first, (*itr).second);
            ++itr;
        }
        printf("\n");
    }

    vector< pair<int, int> > shortestDistances(vertices + 1);
    // shortestDistances is a vector of pairs
    // pair.first -> the shortest distance from start vertex
    // pair.second -> parent vertex in the shortest path

    int startVertex;

    printf("\nEnter a Start Vertex -\n");
    scanf("%d", &startVertex);

    int returnValue = bellmanFord(adjacencyList, vertices, startVertex, shortestDistances);

    if (returnValue == -1)
    {
        printf("No Negative Cycles exist in the Graph -\n");
    }
    else
    {
        printf("Negative Cycles exists in the Graph -\n");
        // The Bellman-Ford Algorithm does not work with negative cycles,
        // all it can do is merely detect them, so when a negative cycle
        // is detected, the shortestDistances vector has wrong values

        PrintNegativeCycle(shortestDistances, shortestDistances[returnValue].second
                           , returnValue);

        return 0;
    }

    printf("\n\nVertex    Shortest Distance to Vertex %d     Parent Vertex-\n", startVertex);
    for (int i = 1; i <= vertices; ++i)
    {
        printf("%d \t  %d \t\t\t\t    %d\n", i, shortestDistances[i].first,
               shortestDistances[i].second);
    }

    return 0;
}




/// Bipartite test

vector <int> adj[mx];
int visited[mx];
int col[mx];
void dfs(int s,int c)
{
    visited[s] = 1;
    col[s] = c ;
    for(auto it : adj[s])
    {

        if(visited[it] == 0)
            {
                if(dfs(it,c^1)==false)
                    return false;
            }
        else
            if(col[s]==col[it])
            return false;
    }
    return true;
}

/// Cycle detection
void dfs(int s,int p)
{
    visited[s] = 1;
    for(auto it : adj[s])
    {
        if(visited[it] == 0)
            {
                if(dfs(it,s)==true)
                    return true;
            }
        else
            if(it==p)
            return true;
    }
    return false;
}


/// Bridges in Graph

vector <int> adj[mx];
bool visited[mx];
int in_time[mx],low_time[mx],timer;

void dfs(int s,int par)
{
    visited[s] = true;
    in_time[s] = low_time[s] = timer;
    timer++;

    for(int i = 0; i < adj[s].size(); ++i)
    {
        if(adj[s][i]==par)
            continue;
        if(visited[adj[s][i]]==true)
        {
            /// back-edge
            low_time[s] = min(low_time[s],in_time[adj[s][i]]);
        }
        else if(visited[adj[s][i]] == false)
            {
                /// forward edge
                dfs(adj[s][i],s);
                if(in_time[s]<low_time[adj[s][i]])
                {
                    /// it is a articulation bridge
                    cout << s << " --> " << adj[s][i] << endl;
                }
                low_time[s] = min(low_time[s],low_time[adj[s][i]]);
            }
    }
}

int main()
{
    int nodes, edges, x, y, connectedComponents = 0;
    cin >> nodes;                       //Number of nodes
    cin >> edges;                       //Number of edges
    for(int i = 0; i < edges; ++i)
    {
        cin >> x >> y;
        //Undirected Graph
        adj[x].push_back(y);                   //Edge from vertex x to vertex y
        adj[y].push_back(x);                   //Edge from vertex y to vertex x
    }

    for(int i = 1; i <= nodes; ++i)
    {
        if(visited[i] == false)
        {
            dfs(i,-1);
        }
    }

}


/// Articulation Point

void Graph::APUtil(int u, bool visited[], int disc[],
                                      int low[], int parent[], bool ap[])
{
    // A static variable is used for simplicity, we can avoid use of static
    // variable by passing a pointer.
    static int time = 0;

    // Count of children in DFS Tree
    int children = 0;

    // Mark the current node as visited
    visited[u] = true;

    // Initialize discovery time and low value
    disc[u] = low[u] = ++time;

    // Go through all vertices aadjacent to this
    list<int>::iterator i;
    for (i = adj[u].begin(); i != adj[u].end(); ++i)
    {
        int v = *i;  // v is current adjacent of u

        // If v is not visited yet, then make it a child of u
        // in DFS tree and recur for it
        if (!visited[v])
        {
            children++;
            parent[v] = u;
            APUtil(v, visited, disc, low, parent, ap);

            // Check if the subtree rooted with v has a connection to
            // one of the ancestors of u
            low[u]  = min(low[u], low[v]);

            // u is an articulation point in following cases

            // (1) u is root of DFS tree and has two or more chilren.
            if (parent[u] == NIL && children > 1)
               ap[u] = true;

            // (2) If u is not root and low value of one of its child is more
            // than discovery value of u.
            if (parent[u] != NIL && low[v] >= disc[u])
               ap[u] = true;
        }

        // Update low value of u for parent function calls.
        else if (v != parent[u])
            low[u]  = min(low[u], disc[v]);
    }
}

// The function to do DFS traversal. It uses recursive function APUtil()
void Graph::AP()
{
    for (int i = 0; i < V; i++)
    {
        parent[i] = NIL;
        visited[i] = false;
        ap[i] = false;
    }

    for (int i = 0; i < V; i++)
        if (visited[i] == false)
            APUtil(i, visited, disc, low, parent, ap);

    /// Now ap[] contains articulation points, print them
    for (int i = 0; i < V; i++)
        if (ap[i] == true)
            cout << i << " ";
}


/// Topological Sort
void dfs(int v) {
    visited[v] = true;
    for (int u : adj[v]) {
        if (!visited[u])
            dfs(u);
    }
    ans.push_back(v);
}

void topological_sort() {
    visited.assign(n, false);
    ans.clear();
    for (int i = 0; i < n; ++i) {
        if (!visited[i])
            dfs(i);
    }
    reverse(ans.begin(), ans.end());
}


/// Warshall
#define V 4
#define INF 99999
void floydWarshall (int graph[][V])
{
	int dist[V][V], i, j, k;

	for (i = 0; i < V; i++)
		for (j = 0; j < V; j++)
			dist[i][j] = graph[i][j];

	for (k = 0; k < V; k++)
		for (i = 0; i < V; i++)
			for (j = 0; j < V; j++)
				if (dist[i][k] + dist[k][j] < dist[i][j])
					dist[i][j] = dist[i][k] + dist[k][j];

}
int main()
{
	int graph[V][V] = { {0, 5, INF, 10},
						{INF, 0, 3, INF},
						{INF, INF, 0, 1},
						{INF, INF, INF, 0}
					};
	floydWarshall(graph);
}



/// MTS - Kruskal

const int MAX = 1e4 + 5;
int id[MAX], nodes, edges;
pair <long long, pair<int, int> > p[MAX];

void initialize()
{
    for(int i = 0;i < MAX;++i)
        id[i] = i;
}

int root(int x)
{
    while(id[x] != x)
    {
        id[x] = id[id[x]];
        x = id[x];
    }
    return x;
}

void union1(int x, int y)
{
    int p = root(x);
    int q = root(y);
    id[p] = id[q];

}

long long kruskal(pair<long long, pair<int, int> > p[])
{
    int x, y;
    long long cost, minimumCost = 0;
    for(int i = 0;i < edges;++i)
    {
        // Selecting edges one by one in increasing order from the beginning
        x = p[i].second.first;
        y = p[i].second.second;
        cost = p[i].first;
        // Check if the selected edge is creating a cycle or not

        if(root(x) != root(y))
        {
            minimumCost += cost;
            union1(x, y);
        }
    }
    return minimumCost;
}

int main()
{
    int x, y;
    long long weight, cost, minimumCost;
    initialize();
    cin >> nodes >> edges;
    for(int i = 0;i < edges;++i)
    {
        cin >> x >> y >> weight;
        p[i] = make_pair(weight, make_pair(x, y));
    }
    // Sort the edges in the ascending order
    sort(p, p + edges);
    minimumCost = kruskal(p);
    cout << minimumCost << endl;
    return 0;
}

/// MST ( prim's )

const int MAX = 1e4 + 5;
typedef pair<long long, int> PII;
bool marked[MAX];
vector <PII> adj[MAX];

long long prim(int x)
{
    priority_queue<PII, vector<PII>, greater<PII> > Q;
    int y;
    long long minimumCost = 0;
    PII p;
    Q.push(make_pair(0, x));
    while(!Q.empty())
    {
        // Select the edge with minimum weight
        p = Q.top();
        Q.pop();
        x = p.second;
        // Checking for cycle
        if(marked[x] == true)
            continue;
        minimumCost += p.first;
        marked[x] = true;
        for(int i = 0;i < adj[x].size();++i)
        {
            y = adj[x][i].second;
            if(marked[y] == false)
                Q.push(adj[x][i]);
        }
    }
    return minimumCost;
}

int main()
{
    int nodes, edges, x, y;
    long long weight, minimumCost;
    cin >> nodes >> edges;
    for(int i = 0;i < edges;++i)
    {
        cin >> x >> y >> weight;
        adj[x].push_back(make_pair(weight, y));
        adj[y].push_back(make_pair(weight, x));
    }
    // Selecting 1 as the starting node
    minimumCost = prim(1);
    cout << minimumCost << endl;
    return 0;
}
