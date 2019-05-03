---
title: '[luogu P3994]高速公路'
date: 2019-04-28 22:59:19
tags:
    - 动态规划
    - 斜率优化
    - 单调栈
categories:
comments:
---

斜率优化dp

<!-- more -->

# 题目背景
C国拥有一张四通八达的高速公路~~网~~树，其中有n个城市，城市之间由一共n-1条高速公路连接。除了首都1号城市，每个城市都有一家本地的客运公司，可以发车前往全国各地，有若干条高速公路连向其他城市，这是一个树型结构，1号城市（首都）为根。假设有一个人要从i号城市坐车出发前往j号城市，那么他要花费Pi*（i城市到j城市的距离）+Qi元。由于距离首都越远，国家的监管就越松，所以距离首都越远，客运公司的Pi（单位距离价格）越大，形式化的说，如果把高速路网看成一棵以首都为根的有根树，i号城市是j号城市的某个祖先，那么一定存在Pi<=Pj。

# 题目描述
大宁成为了国家统计局的调查人员，他需要对现在的高速路网进行一次调查，了解从其他每一个城市到达首都1号城市所花费的金钱(路径必须是简单路径)。

因为有非常多转车（或不转车）的抵达首都的方法，所以人工计算这个结果是十分复杂的。大宁非常的懒，所以请你编写一个程序解决它。
# 题解  
这道题的dp方程比较显然：
$$f[u]=min_{lca(u,v)=v}\{f[v]+P[u]*(depth[u]-depth[v])+Q[u]\}$$
直接dp是$O(n^2)$的。有关于uv的乘积项，depth单调，考虑斜率优化。  
整理得：
$$f[u]=min_{lca(u,v)=v}(-depth[v]*P[u]+f[v])+Q[u]+depth[u]*P[u])$$
则以$-depth[v]$为斜率，$P[u]$为横坐标，$f[v]$为截距，化为一个标准的斜率优化柿子。

我们发现x坐标是不单调的，这意味着我们需要维护整个凸包（单调栈），然后二分找最优的决策点。这也是老套路了。

但是重点并不在于此：这是树上的dp，我们不能像维护序列一样直接令决策点入队出队，因为这样的话每个点不一定只被入队一次，最坏情况仍然是$O(n^2)$。那我们怎么办呢？

~~用主席树实现可持久化栈~~

~~点分治优化dp~~

上面的做法不太好写而且常数巨大...我们有更为优雅的方法：

考虑到斜率单调，那么新来的直线排除掉的旧决策一定是栈顶连续的一段区间。由于决策的单调性，我们可以通过二分找到应该插入新决策的位置(这里找决策点的方法其实和单调队列相同，只是把暴力出队改成二分找合法位置罢了），并让决策入栈。我们发现这样的话其实只是改变了栈顶的位置并修改了一个元素，所以我们在回溯的时候把栈顶和修改的元素改回去，就实现了$O(n\lg n)$的优秀做法...  
~~不过因为数据太水被暴力踩爆了~~

下面是喜闻乐见的代码~（求评价码风qwq）
```cpp
// luogu-judger-enable-o2
#include <cstdio>
#include <algorithm>

const int maxn=1e6+100;

struct Edge
{
    int to,next,w;
}edge[maxn<<1];

int head[maxn],cnt;
int stack[maxn],P[maxn],Q[maxn],fa[maxn];
int64_t f[maxn],depth[maxn];

inline void add(int u,int v,int w)
{
    edge[++cnt].next=head[u];
    edge[cnt].to=v;
    edge[cnt].w=w;
    head[u]=cnt;
}

template<class T>inline T max(T a,T b){return a<b?b:a;}
template<class T>inline T min(T a,T b){return a<b?a:b;}

inline int64_t K(int x){return -depth[x];}
inline int64_t B(int x){return f[x];}
inline int64_t C(int x){return Q[x]+(int64_t)depth[x]*P[x];}
inline double intersection(int x,int y){return ((double)B(x)-B(y))/(K(y)-K(x));}

//f[u]=min_{lca(u,v)=v}(-depth[v]*P[u]+f[v])+Q[u]+depth[u]*P[u])

inline int findbest(int x,int top)
{
    int l=1,r=top;
    while (l<=r)
    {
        int mid=(l+r)>>1;
        if (intersection(stack[mid-1],stack[mid])<=P[x]) l=mid+1;
        else r=mid-1;
    }
    return stack[r];//蒟蒻这里写成了return r WA到怀疑人生..
}

inline int findpos(int x,int top)
{
    int l=1,r=top;
    while (l<=r)
    {
        int mid=(l+r)>>1;
        if (intersection(stack[mid-1],x)>intersection(stack[mid],stack[mid-1])) l=mid+1;
        else r=mid-1;
    }
    return r;
}

inline int64_t calc(int x,int top)
{
    int dec=findbest(x,top);
    return K(dec)*P[x]+B(dec)+C(x);
}

inline void dfs(int u,int top,int ff)
{
    f[u]=calc(u,top);
    top=findpos(u,top)+1;
    int pre=stack[top];
    stack[top]=u;
    for (int i=head[u];i;i=edge[i].next)
    {
        int v=edge[i].to;
        depth[v]=depth[u]+edge[i].w;
        if (v==ff) continue;
        dfs(v,top,u);
    }
    stack[top]=pre;
}

int main()
{
    int n;
    scanf("%d",&n);
    for (int i=2,w;i<=n;++i)
        scanf("%d%d%d%d",fa+i,&w,P+i,Q+i),add(fa[i],i,w);
    for (int i=head[1];i;i=edge[i].next)
    {
        int v=edge[i].to;
        depth[v]=edge[i].w;
        dfs(v,0,1);
    }
    for (int i=2;i<=n;++i)
        printf("%lld\n",f[i]);
}
```