///Linear Diophantine Equation
#include <bits/stdc++.h>
using namespace std;
#define INF 1<<30
#define MAX 10005
#define FASTIO ios_base::sync_with_stdio(false), cin.tie(0), cout.tie(0);
typedef long long ll;

int Extended_gcd(int a, int b, int &x, int &y)
{
    if(a == 0)
    {
        x = 0;
        y = 1;
        return b;
    }
    int x1, y1;
    int d = Extended_gcd(b % a, a, x1, y1);
    x = y1 - (b / a) * x1;
    y = x1;
    return d;
}

//(note that this code does not consider the case a=b=0).We don't consider this case here.
bool find_any_solution(int a, int b, int c, int &x0, int &y0, int &g)
{
    g = Extended_gcd(abs(a),abs(b), x0, y0);
    if(c % g) return false;

    x0 *= c / g;
    y0 *= c / g;
    if(a < 0) x0 = -x0;
    if(b < 0) y0 = -y0;
    return true;
}

void shift_solution(int &x, int &y, int a, int b, int cnt)
{
    x += cnt * b;
    y -= cnt * a;
}

//Note that if a or b is 0, then the problem only has one solution. We don't consider this case here.
int find_all_solution(int a, int b, int c, int minx, int maxx, int miny, int maxy)
{
    int x, y,g;
    if(!find_any_solution(a, b, c, x, y,  g)) return 0;

    a /= g;
    b /= g;

    int sign_a = a > 0 ? +1 : -1;
    int sign_b = b > 0 ? +1 : -1;

    shift_solution(x, y, a, b, (minx - x)/b);
    if(x < minx)
        shift_solution(x, y, a, b, sign_b);
    if(x > maxx)return 0;

    int lx1 = x;

    shift_solution(x, y, a, b, (maxx - x)/b);
    if(maxx > x)
        shift_solution(x, y, a, b, -sign_b);
    int rx1 = x;

    shift_solution(x, y, a,  b, -(miny - y)/a);
    if(y < miny)shift_solution(x, y,  a, b, -sign_a);
    if(y > maxy)
        return 0;
    int lx2 = x;

    shift_solution(x, y, a, b, -(maxy - y)/a);
    if(y > maxy)
        shift_solution(x, y, a, b, sign_a);
    int rx2 = x;

    if(lx2 > rx2) swap(lx2, rx2);

    int lx = max(lx1, lx2);
    int rx = min(rx1, rx2);

    if(lx > rx) return 0;
    return (rx - lx)/ abs(b) + 1;
}


int main()
{

    int T;
    cin >> T;
    for(int cs = 1; cs <= T; cs++){
        int a,b,c,x,y,g;
        cin >> a >> b >> c;
        if(find_any_solution(a,b,c,x,y,g)) cout << "Case "<<cs<<": Yes\n";
        else cout << "Case "<<cs<<": No\n";
    }

    //double end_time = clock();
    //printf( "Time = %lf ms\n", ( (end_time - start_time) / CLOCKS_PER_SEC)*1000);

    return 0;
}
