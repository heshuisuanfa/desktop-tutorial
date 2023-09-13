
memset可以用t代替，要注意t不能为0

#include<iostream>//差分矩阵
#include<cstring>
#include<cstdio>
using namespace std;
int a[1005][1005],d[1006][1005];
void insert(int x1,int y1,int x2,int y2,int x)
{
    d[x1][y1]+=x;
    d[x2+1][y1]-=x;
    d[x1][y2+1]-=x;
    d[x2+1][y2+1]+=x;

}
int main()
{
   int n,m,q;
   cin>>n>>m>>q;
   for(int i=1;i<=n;i++)
   {
       for(int j=1;j<=m;j++)
       {
           cin>>a[i][j];
           insert(i,j,i,j,a[i][j]);
       }
   }
   for(int i=0;i<q;i++)
   {
       int x1,x2,y1,y2,x;
       cin>>x1>>y1>>x2>>y2>>x;
       insert(x1,y1,x2,y2,x);
   }
   for(int i=1;i<=n;i++)
   {
       for(int j=1;j<=m;j++)
       {
           d[i][j]+=d[i-1][j]+d[i][j-1]-d[i-1][j-1];
           printf("%d ",d[i][j]);
       }
       cout<<endl;
   }

    return 0;
}

扩展欧几里得
#include<bits/stdc++.h>
using namespace std;
int exgcd(int a, int b, int &x, int &y){//返回gcd(a,b) 并求出解(引用带回)
    if(b==0){
        x = 1, y = 0;
        return a;
    }
    int x1,y1,gcd;
    gcd = exgcd(b, a%b, x1, y1);
    x = y1, y = x1 - a/b*y1;
    return gcd; 
}
int main(){
    int n,a,b,x,y;
    cin>>n;
    while(n--){
        cin>>a>>b;
        exgcd(a,b,x,y);
        cout<<x<<" "<<y<<endl;
    }
    return 0;
}
1. 扩展欧几里得
用于求解方程 ax+by=gcd(a,b) 的解

当 b=0 时 ax+by=a 故而 x=1,y=0
当 b≠0 时

因为
gcd(a,b)=gcd(b,a%b)

而
bx′+(a%b)y′=gcd(b,a%b)

bx′+(a−⌊a/b⌋∗b)y′=gcd(b,a%b)

ay′+b(x′−⌊a/b⌋∗y′)=gcd(b,a%b)=gcd(a,b)

故而

x=y′,y=x′−⌊a/b⌋∗y′
因此可以采取递归算法 先求出下一层的x′和y′ 再利用上述公式回代即可

2. 对于求解更一般的方程 ax+by=cax+by=c
设 d=gcd(a,b) 则其有解当且仅当 d|c
求解方法如下:

用扩展欧几里得求出 ax0+by0=d 的解

则a(x0∗c/d)+b(y0∗c/d)=c
故而特解为 x′=x0∗c/d,y′=y0∗c/d
而通解 = 特解 + 齐次解

而齐次解即为方程 ax+by=0的解

故而通解为 x=x′+k∗b/d,  y=y′−k∗a/d   k∈z若令 t=b/d，则对于 x 的最小非负整数解为 (x′%t+t)%t
3.应用: 求解一次同余方程 ax≡b(modm)
则等价于求

ax=m∗(−y)+b
ax+my=b


有解条件为 gcd(a,m)|b,然后用扩展欧几里得求解即可

特别的 当 b=1 且 a与m互质时 则所求的x即为a的逆元

线性同余方程
#include<iostream>
#include<cstring>
#include<cstdio>
using namespace std;
int gcd(int a, int b)
{
	return b == 0 ? a : gcd(b, a % b);
}
int n, a, b, m;
int  f(int a, int b, int& x, int& y)
{
	if (b == 0)
	{
		x = 1, y = 0;
		return a;
	}
	int x1, y1, g = gcd(a, b);
	int d = f(b, a % b, x1, y1);
	x = y1, y = x1 - a / b * y1;
	return d;
}
int main()
{
	cin >> n;
	for (int i = 0; i < n; i++)
	{
		cin >> a >> b >> m;
		int g = gcd(a, m);
		if (b % g != 0)
			cout << "impossible" << endl;
		else
		{
			int x, y;
			int d = f(a, m, x, y);
			cout << (long long)x * b/g % m << endl;
		}
	}
	return 0;
}

朴素Dijkstra和堆优化Dijkstra

朴素版dijkstra适合稠密图(边数比较多)
一般用邻接矩阵

寻找路径最短的点：O(n2)O(n2)
加入集合S：O(n)O(n)
更新距离：O(m)O(m)
所以总的时间复杂度为O(n2)

#include<iostream>
#include<cstring>
#include<cstdio>
#include<algorithm>

using namespace std;

const int N = 510, M = 100010;

int g[N][N], dist[N];
bool visited[N];

int n, m;

int dijkstra()
{
    memset(dist, 0x3f, sizeof(dist));
    dist[1] = 0;
    for(int i = 1; i <= n; i++)
    {
        int t = -1;
        for(int j = 1; j <= n; j++)
        {
            if(!visited[j] && (t == -1 || dist[j] < dist[t]))
                t = j;
        }
        visited[t] = true;
        for(int j = 1; j <= n; j++)
            dist[j] = min(dist[j], dist[t] + g[t][j]);
    }
    if(dist[n] == 0x3f3f3f3f) return -1;
    return dist[n];
}

int main()
{
    scanf("%d%d", &n, &m);

    memset(g, 0x3f, sizeof(g));
    while (m--)
    {
        int x, y, c;
        scanf("%d%d%d", &x, &y, &c);
        g[x][y] = min(g[x][y], c);
    }
    cout << dijkstra() << endl;
    return 0;
}

堆优化版本
堆优化版的dijkstra是对朴素版dijkstra进行了优化，在朴素版dijkstra中时间复杂度最高的寻找距离
最短的点O(n^2)可以使用最小堆优化。
1. 一号点的距离初始化为零，其他点初始化成无穷大。
2. 将一号点放入堆中。
3. 不断循环，直到堆空。每一次循环中执行的操作为：
    弹出堆顶（与朴素版diijkstra找到S外距离最短的点相同，并标记该点的最短路径已经确定）。
    用该点更新临界点的距离，若更新成功就加入到堆中。

寻找路径最短的点：O(n)O(n)
加入集合S：O(n)O(n)
更新距离：O(mlogn)

#include<iostream>
#include<cstring>
#include<queue>

using namespace std;

typedef pair<int, int> PII;

const int N = 100010; // 把N改为150010就能ac

// 稀疏图用邻接表来存
int h[N], e[N], ne[N], idx;
int w[N]; // 用来存权重
int dist[N];
bool st[N]; // 如果为true说明这个点的最短路径已经确定

int n, m;

void add(int x, int y, int c)
{
    // 有重边也不要紧，假设1->2有权重为2和3的边，再遍历到点1的时候2号点的距离会更新两次放入堆中
    // 这样堆中会有很多冗余的点，但是在弹出的时候还是会弹出最小值2+x（x为之前确定的最短路径），
    // 并标记st为true，所以下一次弹出3+x会continue不会向下执行。
    w[idx] = c;
    e[idx] = y;
    ne[idx] = h[x]; 
    h[x] = idx++;
}

int dijkstra()
{
    memset(dist, 0x3f, sizeof(dist));
    dist[1] = 0;
    priority_queue<PII, vector<PII>, greater<PII>> heap; // 定义一个小根堆
    // 这里heap中为什么要存pair呢，首先小根堆是根据距离来排的，所以有一个变量要是距离，
    // 其次在从堆中拿出来的时候要知道知道这个点是哪个点，不然怎么更新邻接点呢？所以第二个变量要存点。
    heap.push({ 0, 1 }); // 这个顺序不能倒，pair排序时是先根据first，再根据second，这里显然要根据距离排序
    while(heap.size())
    {
        PII k = heap.top(); // 取不在集合S中距离最短的点
        heap.pop();
        int ver = k.second, distance = k.first;

        if(st[ver]) continue;
        st[ver] = true;

        for(int i = h[ver]; i != -1; i = ne[i])
        {
            int j = e[i]; // i只是个下标，e中在存的是i这个下标对应的点。
            if(dist[j] > distance + w[i])
            {
                dist[j] = distance + w[i];
                heap.push({ dist[j], j });
            }
        }
    }
    if(dist[n] == 0x3f3f3f3f) return -1;
    else return dist[n];
}

int main()
{
    memset(h, -1, sizeof(h));
    scanf("%d%d", &n, &m);

    while (m--)
    {
        int x, y, c;
        scanf("%d%d%d", &x, &y, &c);
        add(x, y, c);
    }

    cout << dijkstra() << endl;

    return 0;
}

bellman_ford
有边数限制的最短路径
#include"iostream"
#include"cstring"
#include"cstdio"
using namespace std;

const int N = 510, M = 10010;

struct Edge {
    int a;
    int b;
    int w;
} e[M];//把每个边保存下来即可
int dist[N];
int back[N];//备份数组防止串联
int n, m, k;//k代表最短路径最多包涵k条边
void bellman_ford() {
    memset(dist, 0x3f, sizeof dist);
    dist[1] = 0;
    for (int i = 0; i < k; i++) {//k次循环
        memcpy(back, dist, sizeof dist);
        for (int j = 0; j < m; j++) {//遍历所有边
            int a = e[j].a, b = e[j].b, w = e[j].w;
            dist[b] = min(dist[b], back[a] + w);
            //使用backup:避免给a更新后立马更新b, 这样b一次性最短路径就多了两条边出来
        }
    }
    if (dist[n] > 0x3f3f3f3f / 2) cout<<"impossible"<<endl;
    else cout<< dist[n]<<endl;

}

int main() {
    scanf("%d%d%d", &n, &m, &k);
    for (int i = 0; i < m; i++) {
        int a, b, w;
        scanf("%d%d%d", &a, &b, &w);
        e[i] = {a, b, w};
    }
    bellman_ford();
    return 0;
}

队列优化bellman_ford

#include<iostream>
#include<cstring>
#include<cstdio>
using namespace std;
const int N = 1e5 + 5;
int n, m, i, j, k;
int u[N], v[N], w[N];
int first[N], ne[N], dis[N], book[N], que[N], head = 1, tail = 1;
int inf = 999999999;

int main()
{
	scanf("%d%d", &n, &m);
	for (i = 1; i <= n; i++)
		dis[i] = inf;

	dis[1] = 0;

	for (i = 1; i <= n; i++)
		first[i] = -1;

	for (i = 1; i <= m; i++)
	{
		scanf("%d%d%d", &u[i], &v[i], &w[i]);
		ne[i] = first[u[i]];
		first[u[i]] = i;
	}
	que[tail++] = 1;
	book[1] = 1;
	while (head < tail)
	{
		k = first[que[head]];
		while (k != -1)
		{
			if (dis[v[k]] > dis[u[k]] + w[k])
			{
				dis[v[k]] = dis[u[k]] + w[k];
				if (book[v[k]] == 0)
				{
					que[tail++] = v[k];
					book[v[k]] = 1;
				}
			}
			k = ne[k];
		}
		book[que[head]] = 0;
		head++;
	}
	for (i = 1; i <= n; i++)
	{
		printf("%d ", dis[i]);
	}
	return 0;
}

**数据构造**

```c++
#include<iostream>
#include<time.h>
#include<stdlib.h>
using namespace std;
int main()
{
	srand(time(nullptr));
	int t = 100;
	while (t--)
	{
		int n = rand() % 10 + 1;//产生一个随机数
		cout << n << endl;
		for (int i = 0; i < n; i++)
		{
			cout << rand() % 1000 << ' ';
		}
		cout << endl;
		for (int i = 0; i < n; i++)
		{
			cout << rand() % 2 << ' ';
		}
		cout << endl;
	}

	return 0;
}
```
