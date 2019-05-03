---
title: '[SCOI2014][Luogu P3288]方伯伯运椰子'
date: 2019-04-22 22:44:36
tags:
categories:
comments:
---
# 题目描述
四川的方伯伯为了致富，决定引进海南的椰子树。方伯伯的椰子园十分现代化，椰子园中有一套独特的交通系统。

<!-- more -->

现在用点来表示交通节点，边来表示道路。这样，方伯伯的椰子园就可以看作一个有 n + 2 个交通节点，m条边的有向无环图。n +1 号点为入口，n +2 号点为出口。每条道路都有 6 个参数，ui，vi，ai，bi，ci，di，分别表示，该道路从 ui 号点通向 vi 号点，将它的容量压缩一次要 ai 的花费，容量扩大一次要 bi 的花费，该条道路当前的运输容量上限为 ci，并且每单位运输量通过该道路要 di 的费用。

在这个交通网络中，只有一条道路与起点相连。因为弄坏了这条道路就会导致整个交通网络瘫痪，聪明的方伯伯决定绝不对这条道路进行调整，也就是说，现在除了这条道路之外，对其余道路都可以进行调整。

有两种调整方式：

选择一条道路，将其进行一次压缩，这条道路的容量会下降 1 单位。

选择一条道路，将其进行一次扩容，这条道路的容量会上升 1 单位。

一条道路可以被多次调整。

由于很久以前，方伯伯就请过一个工程师，对这个交通网络进行过一次大的优化调整。所以现在所有的道路都被完全的利用起来了，即每条道路的负荷都是满的（每条道路的流量等于其容量）。

但方伯伯一想到自己的海南椰子会大丰收，就十分担心巨大的运输量下，会导致过多的花费。因此，方伯伯决定至少进行一次调整，调整之后，必须要保持每条道路满负荷，且总交通量不会减少。

设调整后的总费用是 Y，调整之前的总费用是 X。现在方伯伯想知道，最优调整比率是多少，即假设他进行了 k 次调整，(X - Y)/k最大能是多少？

注：总费用 = 交通网络的运输花费 + 调整的花费

# 输入输出格式
## 输入格式：
第一行包含二个整数N，M接下来M行代表M条边，表示这个交通网络每行六个整数，表示Ui,Vi,Ai,Bi,Ci,Di接下来一行包含一条边，表示连接起点的边

## 输出格式：
一个浮点数，保留两位小数。表示答案，数据保证答案大于0  
# 数据规模与约定
```
1<=N<=5000
0<=M<=3000
1<=Ui,Vi<=N+2
0<=Ai,Bi<=500
0<=Ci<=10000
0<=Di<=1000
```

# 题解
~~真是个语文题~~  
在艰难看懂了题意之后，我们可以发现，方伯伯实际上是在DAG上跑网络流。所以问题可以简化为：给一个网络图，可以进行增广或者退流，对边i增广一次的费用为$b_i+d_i$,退流一次的费用为$a_i-d_i$,（注意边权为0的边不能退流）最终使得总流量不变（增广和退流都需要花费，所以总流量变大肯定不优），并最大化原式。

题面给的$(X-Y)/k$太难看了，我们给它换一下：
令$Y=X+\Delta w$，
则原式变为$-\Delta w/k$。

如何保持流量不变的情况下，优化增广的费用？

我们需要一个定理：
> 消圈定理：一个流是当前流量下的最小费用流，等价于当前残量网络上没有负费用圈。

> 证明：假设残量网络上存在负费用圈，我们可以把原来不经过负圈的流沿着这个负圈增广一次，则我们的流量会不变，并且减少了所花的费用。

所以：我们现在的目标就是在残量网络上沿着负圈增广,且$\text{maximize}\{-\Delta w/k\}$.
因为我们要最大化比值，所以我们显然要选择**一个**最优的环进行增广。选择多个环，答案显然会变劣（求的是平均值嘛）。这就是一个最小（大）平均值回路问题了。套用分数规划的思路，二分+bellman-ford判断负环就可以了。
# 代码
```cpp
// luogu-judger-enable-o2
#include <cstdio>
#include <cstring>
#include <queue>

using std::queue;

const int maxn=5000+1000;
const double eps=1e-6;

struct Edge
{
    int to,next;double w;
}edge[maxn<<1];

int head[maxn],cnt;
int u[maxn],v[maxn],a[maxn],b[maxn],c[maxn],d[maxn];

inline void add(int u,int v,double w)
{
    edge[++cnt].next=head[u];
    edge[cnt].to=v;
    edge[cnt].w=w;
    head[u]=cnt;
}

bool spfa(int s,int n)
{
    static bool vis[maxn];
    static double d[maxn];
    static int cnt[maxn];
    queue<int> q;
    q.push(s);
    for (int i=1;i<=n;++i) d[i]=1e18;
    memset(vis,0,sizeof(vis));
    memset(cnt,0,sizeof(cnt));
    vis[s]=true;d[s]=0;
    while (!q.empty())
    {
        int u=q.front();q.pop();
        vis[u]=false;
        for (int i=head[u];i;i=edge[i].next)
        {
            int v=edge[i].to;
            if (d[v]>d[u]+edge[i].w)
            {
                d[v]=d[u]+edge[i].w;
                cnt[v]=cnt[u]+1;
                if (cnt[v]>=n) return true;
                if (!vis[v]) vis[v]=true,q.push(v);
            }
        }
    }
    return false;
}

inline bool check(double mid,int n,int m)
{
    cnt=0;
    memset(head,0,sizeof(head));
    for (int i=1;i<=m;++i)
    {
        add(u[i],v[i],mid+(b[i]+d[i]));
        if (c[i]>0) add(v[i],u[i],mid+(a[i]-d[i]));
    }
    return spfa(n-1,n);
}

int main()
{
    int n,m;
    scanf("%d%d",&n,&m);
    n+=2;
    double sigma=0;
    for (int i=1;i<=m;++i)
        scanf("%d%d%d%d%d%d",u+i,v+i,a+i,b+i,c+i,d+i),sigma+=a[i]+b[i]+(double)c[i]*d[i];
    double l=0,r=sigma;
    while (l+eps<r)
    {
        // printf("%.3lf %.3lf\n",l,r);
        double mid=(l+r)/2;
        if (check(mid,n,m)) l=mid;
        else r=mid;
    }
    printf("%.2lf",l);
}
```