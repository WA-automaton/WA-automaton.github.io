---
title: YNOI做题记录
date: 2019-02-06 21:43:29
tags:
    - 分块
    - 莫队
    - YNOI
    - 毒瘤
    - 数据结构
    - 平衡树
categories:
    - 题解
comments:
---

> 万古神犇LXL，数据结构碾众生！

~~即使是蒟蒻也想变强啊..~~

# [luogu P3987](https://www.luogu.org/problemnew/show/P3987) 我永远喜欢珂朵莉~

## 题目大意：

有一个长为n的非负数序列A，支持以下两个操作：

+ 1 l r x : 把区间[l,r]中所有x的倍数/x
+ 2 l r : 查询区间[l,r]的和

## 数据范围：

$1 \le n , m \le 100000$

$0 \le A_i \le 500000$

$1 \le x \le 500000$

<!-- more-->

## 题解：

首先，这道题的突破口在这里：

* 一个数的约数个数不会太多。虽说上界$O(\sqrt n)$,但实际上远没有那么多。500000以内只有大概200个左右。当对值域进行限制的时候，很可能就与约数个数相关。尤其是，本题中还有很明显的x的倍数÷x的操作，珂以考虑对每个可能的约数进行维护。

  鉴于单点修改，区间求和是$O(\log n)$ 的，而一个数最多被除$O(\log n)$次，总复杂度$O(n\log^2 n)$这并不是制约复杂度的关键。而且，如果对于整个数列或分块以后的整块（本质上是一个“整体”）维护整体的信息的话，这题根本不可做了。必须找出需要被除的数，才能维护整个数列。现在问题是，如何快速找到需要被除的数。即，如何快速找到序列中x的倍数。珂以考虑对所有$x\in (2,500000)$进行维护。对于每个约数，维护一棵平衡树，存储数列中它的倍数的下标。当进行区间除的时候，在对应的平衡树中找到下标在$[l, r]$之间的子树，进行dfs，并删掉所有操作进行后不再是x倍数的数的下标。这个过程中珂以顺便维护数列值的变化。当然，并不需要建出所有的平衡树，只对查询的x建树就行了。

  说来惭愧，我这道题在看了lxl的题解后还改了好几天。我从这道题吸取的经验有以下几点：

  + 当你的板子检查了很多很多遍都没发现问题时，很可能是main函数写错了（捂脸
  + 各种最大值一定要弄清楚，例如我写题的时候就把值域最大值当成了n。
  + 板子是珂以根据自己的需要而改动的。如本题中平衡树并不需要维护size。

当然YNOI的毒瘤题需要一点卡常的小trick，相信大家都会，不再赘述。

### 我的代码：

```cpp
#include <cmath>
#include <cstdio>
#include <cctype>
#include <vector>
#include <algorithm>
#include <cassert>
#include <climits>
#include <ctime>
#define ls(o) (t[o].ch[0])
#define rs(o) (t[o].ch[1])

using std::vector;

typedef long long ll;

const int maxn=500000+1000;
const int INF=0x3f3f3f3f;

int a[maxn/5],cnt,n,root[maxn],cntdel,del[maxn];
vector<int> v[maxn];
ll c[maxn];

namespace IO
{
    static char buf[1<<25],*fs,*ft;
    // inline char gc()
    // {
    //     if (fs==ft)
    //     {
    //         ft=(fs=buf)+fread(buf,1,1<<25,stdin);
    //         if (fs==ft) return EOF;
    //     } 
    //     return *fs++;
    // }
    // #define gc() getchar()
    // #define gc (*fs++)
    inline int read()
    {
        register char ch;
        while (!isdigit(ch=(*fs++)));
        register int x=ch-48;
        while (isdigit(ch=(*fs++)))
            x=x*10+ch-48;
        return x;
    }
}using IO::read;

struct Node
{
    int ff,ch[2];
    int pos;
}t[maxn*70];

template<class T> inline T max(T a,T b){return a<b?b:a;}
template<class T> inline T min(T a,T b){return a<b?a:b;}

inline void update(int p,int x)
{
    for (register int i=p;i<=n;i+=i&-i) c[i]+=x;
}

inline ll query(int p)
{
    ll ans=0;
    for (register int i=p;i;i-=i&-i) ans+=c[i];
    return ans;
}

inline void decompose(int x,int p)
{
    int n=sqrt(x)+0.1;
    for (register int i=1;i<=n;++i)
    {
        if (x%i==0)
        {
            v[i].push_back(p);
            if (i*i!=x) v[x/i].push_back(p);
            // printf("%d %d",i,x/i);
        }
    }
}

inline void rotate(int x)
{
    int y=t[x].ff,z=t[y].ff;
    int k=t[y].ch[1]==x;
    t[z].ch[t[z].ch[1]==y]=x;
    t[x].ff=z;
    t[y].ch[k]=t[x].ch[k^1];
    t[t[x].ch[k^1]].ff=y;
    t[x].ch[k^1]=y;
    t[y].ff=x;
}

int build(int l,int r,const vector<int>& v,int fa)
{
    if (l>r) return 0;
    int mid=(l+r)>>1;
    int u=++cnt;
    t[u].ch[0]=build(l,mid-1,v,u);
    t[u].ch[1]=build(mid+1,r,v,u);
    t[u].ff=fa;
    t[u].pos=v[mid];
    return u;
}

int _pre,_suc,ppre,psuc;

void _prec(int x,int o)
{
    if (!o) return;
    if (t[o].pos<x && _pre<t[o].pos)
        _pre=t[o].pos,ppre=o;
    _prec(x,t[o].ch[t[o].pos<x]);
}

inline int prec(int x,int o)
{
    _pre=INT_MIN;ppre=0;
    _prec(x,o);
    return ppre;
}

void _succ(int x,int o)
{
    if (!o) return;
    if (t[o].pos>x && _suc>t[o].pos)
        _suc=t[o].pos,psuc=o;
    _succ(x,t[o].ch[t[o].pos<=x]);
}

inline int succ(int x,int o)
{
    _suc=INT_MAX;psuc=0;
    _succ(x,o);
    return psuc;
}

inline void splay(int x,int goal,int idx)
{
    while (t[x].ff!=goal)
    {
        int y=t[x].ff,z=t[y].ff;
        if (z!=goal) rotate((t[z].ch[1]==y)^(t[y].ch[1]==x)?x:y);
        rotate(x);
    }
    if (!goal) root[idx]=x;
}

void remove(int x,int u,int idx)
{
    // assert(t[x].pos>-INF && t[x].pos<INF);
    int L=prec(t[x].pos,u),R=succ(t[x].pos,u);
    if (L!=R)
    {
        splay(L,0,idx);
        splay(R,root[idx],idx);
        // assert(t[t[t[root[idx]].ch[1]].ch[0]].pos==t[x].pos && !t[x].ch[0] && !t[x].ch[1]);
        t[t[t[root[idx]].ch[1]].ch[0]].ff=0;
        t[t[root[idx]].ch[1]].ch[0]=0;
    }
}

void dfs(int u,int x,int idx)
{
    if (!u) return;
    if (ls(u)) dfs(ls(u),x,idx);
    if (rs(u)) dfs(rs(u),x,idx);
    // assert(a[t[u].pos]%x==0);
    if (a[t[u].pos]%x!=0) del[++cntdel]=u;
    else
    {
        update(t[u].pos,-a[t[u].pos]);
        a[t[u].pos]/=x;
        update(t[u].pos,a[t[u].pos]);
        if (a[t[u].pos]%x!=0) del[++cntdel]=u;
    }
}

void Divide(int l,int r,int x)
{
    int L=prec(l,root[x]),R=succ(r,root[x]);
    assert(R && L);
    // assert(root[x]);
    if (L!=R)
    {
        splay(L,0,x);
        splay(R,root[x],x);
        cntdel=0;
        dfs(t[t[root[x]].ch[1]].ch[0],x,x);
        // for (int i=1;i<=cntdel;++i)
        //     assert(t[del[i]].pos>=l && t[del[i]].pos<=r);
        for (int i=1;i<=cntdel;++i) remove(del[i],root[x],x);
    }
}

void printtree(int u)
{
    if (ls(u)) printtree(ls(u));
    printf("%d ",t[u].pos);
    // if (t[u].pos>-INF && t[u].pos<INF) assert(a[t[u].pos]%x==0);
    if (rs(u)) printtree(rs(u));
}

int main()
{
    // freopen("input.in","r",stdin);
    // freopen("my.out","w",stdout);
    fread(IO::fs=IO::buf,1,1<<25,stdin);
    register int m;
    n=read();m=read();
    for (register int i=2;i<=500000;++i) v[i].push_back(-INF);
    for (register int i=1;i<=n;++i)
        a[i]=read(),decompose(a[i],i),update(i,a[i]);
    // for (register int i=2;i<=500000;++i)
    //     if (v[i].size()>1) v[i].push_back(+INF),root[i]=build(0,v[i].size()-1,v[i],0);
    // for (int i=2;i<=n;++i)
    //     if (v[i].size()>2)printtree(root[i],i);
    static bool used[maxn];
    for (register int i=1,opt,l,r,x;i<=m;++i)
    {
        opt=read();l=read();r=read();
        if (opt==1)
        {
            x=read();
            if (!used[x] && v[x].size()>1 && x>1)
            {
                v[x].push_back(+INF);
                root[x]=build(0,v[x].size()-1,v[x],0);
                used[x]=true;
            }
            if (l>r) l^=r^=l^=r;
            if (l<0 || r>n) continue;
            if (x>1 && v[x].size()>2) Divide(l,r,x);
        }
        else printf("%lld\n",query(r)-query(l-1));
    }
    // freopen("tm.out","w",stdout);
    // printf("%d",clock());
    // for(;;);
}
```

# [luogu P5068](https://www.luogu.org/problemnew/show/P5068)[Ynoi2015]我回来了

## 题目描述

珂朵莉给你一个无向图，每次查询的时候给一堆二元组$(x_i,y_i)$
求图中有多少个点u与至少一个这次询问给出的二元组$(x_i,y_i)$满足 $dist(u,x_i)\le y_i$，dist表示这两个点在图中的距离
如果不连通$dist = \inf$

## 输入输出格式

### 输入格式：

* 第一行三个数表示n，m，q
* n表示顶点个数，m表示边数
* 之后m行每行两个数x，y表示这两个点之间连有一条边~，边权都为1
* 之后q次询问，每个询问先给你一个数a
* 之后a行每行两个数，x，y，表示一个二元组
* n <= 1000 , m <= 100000 , q <= 100000
* a的和 <= 2100000

### 输出格式：

* q行，每行一个数表示这次询问的答案

## 题解：

这道题大概是YNOI中最良心的一道题？(雾  ~~然而我也没做出来~~

首先，我们珂以想到，这题大概要预处理，然后用接近$O(1)$的时间回答每个二元组询问。考虑题目中，每一个点即使满足所有要求也只被计算一次贡献，而每个二元组之间又是相互独立的，所以我们需要一个资瓷快速集合取并，快速求集合元素数目的数据结构——这不就是$\texttt{bitset}$嘛。所以，我们尝试，令$f(u,i)$为$dist(u,v)\le i$的v的集合，则对于每次询问，将所对应的$f$集合取一个并就好了。现在，问题转化为如何求f集合。本题的时空限制十分诡异，会让人误以为标算的复杂度是错的~~（然鹅实际跑得飞快）~~其实暴力计算f就珂以了~~

* 首先，跑n次bfs，计算出所有顶点对之间最短路。时间$O(n(n+m))$
* 然后，用已知的信息暴力更新f。设w为bitset的位宽，时间$O(n^3/w)$
* （逃

最后，毒瘤lxl卡了链式前向星，只能用vector存图（因为空间够用，所以我直接开了数组）

注意，输入存在重边。

所以，这道题我们就做完啦（撒花！

### 我的代码

~~（我以为这题一定很卡常，所以代码稍微毒瘤了一点）~~：

```cpp
#include <cstdio>
#include <cctype>
#include <cstring>
#include <bitset>
#include <queue>

using std::queue;
using std::bitset;

const int maxn=1024;
const int INF=0x3f3f3f3f;

int G[maxn][maxn],cnt[maxn],d[maxn][maxn];
int link[maxn][maxn];

namespace IO
{
    static char buf[1<<17],*fs,*ft;
    inline char gc()
    {
        if (fs==ft)
        {
            ft=(fs=buf)+fread(buf,1,1<<17,stdin);
            if (fs==ft) return EOF;
        }
        return *fs++;
    }
    #define gc getchar
    inline int read()
    {
        char ch;
        while (!isdigit(ch=gc()));
        int x=ch^48;
        while (isdigit(ch=gc()))
            x=x*10+ch-48;
        return x;
    }
}-using IO::read;

inline void bfs(int s)
{
    static bool vis[maxn];
    memset(vis,0,sizeof(vis));
    memset(d[s],0x3f,sizeof(int)*maxn);
    vis[s]=true;d[s][s]=0;
    queue<int> q;
    q.push(s);
    while (!q.empty())
    {
        register int u=q.front();q.pop();
        for (register int *ptr=G[u],*ed=G[u]+cnt[u];ptr!=ed;++ptr)
        {
            register int v=*ptr;
            if (!vis[v]) d[s][v]=d[s][u]+1,q.push(v),vis[v]=true;
        }
    }
}

int main()
{
    int n,m,q;
    n=read();m=read();q=read();
    for (register int i=1,u,v;i<=m;++i)
    {
        u=read();v=read();
        if (link[u][v]) continue;
        G[u][cnt[u]++]=v;
        G[v][cnt[v]++]=u;
        link[u][v]=link[v][u]=true;
    }
    for (register int i=1;i<=n;++i) bfs(i);
    static bitset<maxn> f[maxn][maxn];
    for (register int u=1;u<=n;++u)
        for (register int v=1;v<=n;++v)
            if (d[u][v]<=n) f[u][d[u][v]].set(v);
    for (register int u=1;u<=n;++u)
        for (register int i=1;i<=n;++i)
            f[u][i]|=f[u][i-1];
    bitset<maxn> ans;
    for (register int i=1;i<=q;++i)
    {
        ans.reset();
        register int a=read();
        for (register int j=1;j<=a;++j)
        {
            register int x=read(),y=read();
            ans|=f[x][y];
        }
        printf("%d\n",(int)ans.count());
    }
}
```



# [luogu P5072](https://www.luogu.org/problemnew/show/P5072)[YNOI2015]盼君勿忘

## 题目描述

珂朵莉给了你一个序列，每次查询一个区间[l,r]中所有子序列分别去重后的和mod p

## 输入输出格式

### 输入格式：

* 第一行两个数$n,m$
* 第二行$n$个数表示这个序列$A$
* 之后$m$行，每行三个数$l,r,p$表示查询的区间与模数
* 对于$100\%$的数据，$n,m <= 10^5$,$A_i,p <= 1000000000$

### 输出格式：

* $m$行，每行输出一个数表示答案

## 题解

Orz lxl...由乃OI真的毒瘤。。。

~~我从不抄题解，我只是题解的搬运工。~~真的不会做，瞎写一点儿吧。

看到这题，没修改操作，数据范围1e5，珂以想到这是个莫队题。~~然后我就不会了~~。

看到这样的题，我们暴力考虑每个区间显然是不现实的。我们转而考虑每个数对答案的贡献。对于一个数，它被计算到的次数珂以这么计算：

对于一个区间$[l..r]$，它的非空子序列个数显然为$$2^{r-l+1}-1$$个。不包含x的子序列有$$2^{r-l+1-cnt[x]}-1$$个。所以，x一共被计算了$2^{r-l+1}-2^{r-l+1-cnt[x]}$次。

~~所以，我们用$cnt[x]$维护x在当前区间的出现次数，莫队维护cnt，修改的同时暴力维护和。。。这就珂以做到$O(n\sqrt n)$了~~

上面的假做法只适合模数不变的情况。。。我们考虑模数会变化该怎么做。

**手动划重点！以下内容的思想很重要！**

考虑以$\sqrt n$为界限**分类讨论**。对于出现次数小于$\sqrt n$的数，他们的出现次数只有$\sqrt n$种。对于出现次数大于等于$\sqrt n$的数，这样的数最多只有$\sqrt n$种。所以我们对于每个小于等于$\sqrt n$ 的出现次数x，维护出现次数=x的数之和，对于出现次数$\geq \sqrt n$的数，用一种数据结构维护它们组成的集合。$\texttt{STL unordeded_set}$是个好东西啊。~~虽然考试不让用~~这样就珂以做到$O(n\sqrt n)$求解问题啦。

最后还有一个小问题——如何快速求2的幂？模数动态改变，每次必须重新计算。有一种$O(\sqrt n)$预处理，$O(1)$查询的处理方法，正好适合本题。

设$T=\lfloor\sqrt n\rfloor$,则$C^{x}=C^{\lfloor x/T \rfloor} * C^{x\mod T}$ 。这样，只需要预处理$C^1, C^2, C^3....C^T$和$C^T, C^{2T}....C^{T*T}$ 即可。

### 我的代码：

(卡常的血泪史。。。)另外，好像对于值域比较大的情况，它会输出负数，现在暂时没找到原因。还有，膜运算真的太慢啦！去掉四个膜运算，运行时间几乎-1s。。

```cpp
#pragma GCC optimize("-O3")
#pragma GCC optimize("-Ofast")
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <ctime>
// #define int long long

typedef unsigned long long ll;

using std::sort;
using std::unique;
using std::lower_bound;
using std::unordered_map;
using std::unordered_set;

const int maxn=1e6+100;

int blo[maxn],a[maxn],b[maxn],cnt[maxn],L=1,R;
ll _2_k[maxn],_2_sqrt[maxn];
ll sum[maxn],ans;
unordered_map<int,int> rnk;
unordered_set<int> mp;

namespace IO
{
    char buf[1<<24],*fs,*ft;
    // inline char gc()
    // {
    //     if (fs==ft)
    //     {
    //         ft=(fs=buf)+fread(buf,1,1<<24,stdin);
    //         if (fs==ft) return EOF;
    //     }
    //     return *fs++;
    // }
    #define gc() (*fs++)
    // #define gc() getchar()
    inline int read()
    {
        char ch;
        while (!isdigit(ch=gc()));
        int x=ch^48;
        while (isdigit(ch=gc()))
            x=x*10+ch-48;
        return x;
    }
}using IO::read;

struct Q
{
    int l,r,p,id;
    bool operator< (const Q& b) const 
	{
		// return blo[l]^blo[b.l]?blo[l]<blo[b.l]:blo[l]&1?r<b.r:r>b.r;
        return blo[l]==blo[b.l]?r<b.r:l<b.l;
	}
}q[maxn];

inline ll qpow(register ll a,register ll b,register ll p)
{
    ll ans=1%p;
    for (;b;b>>=1)
    {
        if (b&1) ans=ans*a%p;
        a=a*a%p;
    }
    return ans;
    // fprintf(stderr,"SYTAKIOI\n");
}

inline void ins(const int t,const int T)
{
    register int &x=cnt[t];
    if (x<T)
    {
        if (x) sum[x]-=b[t];
        if (++x<T) sum[x]+=b[t];
        else mp.insert(b[t]);
    }
    else ++cnt[t];
    // fprintf(stderr,"SYTAKIOI\n");
}

inline void del(const int t,const int T)
{
    register int &x=cnt[t];
    register int y=b[t];
    if (x<T)
    {
        sum[x]-=y;
        if (--x) sum[x]+=y;
    }
    else if (--x<T) mp.erase(y),sum[T-1]+=y;
            // fprintf(stderr,"SYTAKIOI\n");
}

inline void prework(const ll p,const int T)
{
    _2_k[0]=1%p,_2_k[1]=2%p,_2_sqrt[0]=1%p,_2_sqrt[1]=qpow(2,T,p);
    for (int i=2;i<=T;++i)
    {
        _2_k[i]=_2_k[i-1]*2%p;
        _2_sqrt[i]=_2_sqrt[i-1]*_2_sqrt[1]%p;
    }
    // fprintf(stderr,"SYTAKIOI\n");
}

inline ll getpow(const int x,const int T,const ll p)
{
    // fprintf(stderr,"SYTAKIOI\n");
    return _2_k[x%T]*_2_sqrt[x/T]%p;
}

signed main()
{
    freopen("in.txt","r",stdin);
    freopen("my.out","w",stdout);
    fread(IO::fs=IO::buf,1,1<<24,stdin);
    register int n,m;
    n=read();m=read();
    const int T=sqrt(n);
    for (register int i=1;i<=n;++i)
        blo[i]=(i-1)/T+1;
    for (register int *ptrb=b+1,*ptra=a+1,*ed=a+n+1;ptra!=ed;++ptrb,++ptra)
        *ptrb=*ptra=read();
    sort(b+1,b+n+1);
    int tot=unique(b+1,b+n+1)-b-1;
    for (register int i=1;i<=n;++i)
        a[i]=lower_bound(b+1,b+tot+1,a[i])-b;
    // for (int i=1;i<=n;++i)
    //     printf("x=%d rnk[x]=%d\n",a[i],rnk[a[i]]);
    for (register Q* ptr=q+1,*ed=q+m+1;ptr!=ed;++ptr)
        ptr->l=read(),ptr->r=read(),ptr->p=read(),ptr->id=ptr-q;
    sort(q+1,q+m+1);
    static ll ans[maxn];
    for (register Q* ptr=q+1,*ed=q+m+1;ptr!=ed;++ptr)
    {
        // if (q[i].p==15) q[i].p=1e9+7;
        const int l=ptr->l,r=ptr->r,p=ptr->p;
        prework(p,T);
        while (R<r) ins(a[++R],T);
        while (L>l) ins(a[--L],T);
        while (R>r) del(a[R--],T);
        while (L<l) del(a[L++],T);
        ll& tot=ans[ptr->id];
        ll tmp=getpow(r-l+1,T,p);
        for (register int i=1;i<T;++i)
            (tot+=((tmp-getpow(r-l+1-i,T,p)%p)+p)*sum[i])%=p;
        for (const auto x:mp)
            (tot+=((tmp-getpow(r-l+1-(cnt[x]),T,p)%p)+p)*b[x])%=p;
        // printf("L=%d R=%d,sum[1]=%d\n",L,R,sum[1]);
    }
    for (register ll *i=ans+1,*ed=ans+m+1;i!=ed;++i)
        printf("%llu\n",*i);
    // fprintf(stderr,"%d",clock());
    // return clock();
}
```



# [luogu P3674](https://www.luogu.org/problemnew/show/P3674)小清新人渣的本愿

## 题目描述

给你一个序列a，长度为n，有m次操作，每次询问一个区间是否珂以选出两个数它们的差为x，或者询问一个区间是否珂以选出两个数它们的和为x，或者询问一个区间是否珂以选出两个数它们的乘积为x ，这三个操作分别为操作1,2,3

选出的这两个数珂以是同一个位置的数

## 输入输出格式

### 输入格式：

* 第一行两个数$n,m$
* 后面一行$n$个数表示$a_i$
* 后面$m$行每行四个数$opt\ l\ r\ x$
* $opt$表示这个是第几种操作，$l,r$表示操作的区间，$x$表示这次操作的$x$

### 输出格式：

* 对于每个询问，如果珂以，输出$\texttt{hana}$，否则输出$\texttt{bi}$

## 题解

这题还是挺良心的~~比上一题良心多了~~  ~~瞎说毒瘤lxl怎么会出良心题~~

没修改+1e5+lxl=膜队。

所以我们考虑怎么膜队维护。

我们首先肯定要维护所有数的出现次数。然后怎么做？

考虑暴力：对于每个可能成为答案的数，查询对应的另一个数是否存在。这个大概没法用什么数据结构优化了，所以考虑神奇的$\texttt{bitset}$。用$\texttt{bitset}$维护每个数是否出现，我们发现：

* 对于询问1，若存在k和k+x，则对应着$\texttt{S&(S<<x).any()==true}$ 。
* 对于询问2，若存在$a+b=x$,即$a-(-b)=x$。对于负数，我们不好维护，所以用处理负下标的一般方法，转化成$a-(-b+N)=x-N$。（此处N是一个较大的整数）这样，用另外一个$\texttt{bitset}$，对于每个出现的$k$，令对应的$(-k+n)$下标处的值为1就珂以了。另外，$\texttt{bitset}$的移位运算会将$\texttt{int}$强转成$\texttt{size_t}$，当移位数为负数的时候会出锅。所以要手动把上面的左移改成右移，移位数取反。
* 对于询问3，我们发现我们用$\texttt{bitset}$没法很好的维护了。所以暴力枚举因子，判断是否存在,珂以发现这并不是复杂度的瓶颈。这样做是正确的......

### 我的代码：

~~您看我如果不卡常的话码风还是挺正常的嘛~~

```cpp
#include <cmath>
#include <bitset>
#include <cstdio>
#include <algorithm>

using std::bitset;
using std::sort;

const int maxn=1e5+1000;

int blo[maxn],a[maxn],col[maxn];
bitset<maxn> ex,rev;
int d[maxn];

struct Q
{
    int type,l,r,x,id;
    bool operator< (const Q& q) const
    {
        return blo[l]==blo[q.l]?r<q.r:l<q.l;
    }
}q[maxn];

inline bool QuerySub(int x)
{
    return (ex&(ex<<x)).any();
}

inline bool QueryAdd(int x)
{
    return (ex&(rev>>(maxn-x))).any();
}

inline bool QueryMul(int x)
{
    for (int i=1;i*i<=x;++i)
    {
        if (x%i==0)
            if (ex[i] && ex[x/i]) return true;
    }
    return false;
}

inline void add(int x)
{
    if (++col[x]==1) 
        ex[x]=1,rev[maxn-x]=1;
}

inline void del(int x)
{
    if (--col[x]==0)
        ex[x]=0,rev[maxn-x]=0;
}

int main()
{
    int n,m;
    scanf("%d%d",&n,&m);
    int T=sqrt(n);
    for (int i=1;i<=n;++i)
        scanf("%d",a+i),blo[i]=(i-1)/T+1;
    for (int i=1;i<=m;++i)
    {
        scanf("%d%d%d%d",&q[i].type,&q[i].l,&q[i].r,&q[i].x);
        q[i].id=i;
    }
    sort(q+1,q+m+1);
    static int ans[maxn];
    for (int i=1,L=1,R=0;i<=m;++i)
    {
        const int l=q[i].l,r=q[i].r,x=q[i].x;
        while (R<r) add(a[++R]);
        while (L>l) add(a[--L]);
        while (L<l) del(a[L++]);
        while (R>r) del(a[R--]);
        switch (q[i].type)
        {
            case 1:
                ans[q[i].id]=QuerySub(x);
                break;
            case 2:
                ans[q[i].id]=QueryAdd(x);
                break;
            case 3:
                ans[q[i].id]=QueryMul(x);
                break;
        }
    }
    for (int i=1;i<=m;++i) puts(ans[i]?"hana":"bi");
}
```



# [luogu P4688](https://www.luogu.org/problemnew/show/P4688)[YNOI2016]掉进兔子洞

## 题目描述

给定一个长度为n的数列，有m次询问。每次询问给出三个区间，询问若从这三个区间一起删数，一直删到三个区间没有共同的数为止，最后这三个区间一共剩下的数的数量。（询问之间互相独立）

## 题解

这道题的话，首先，还是经验公式：

没修改+1e5+lxl=膜队。

所以考虑怎么膜队。

首先通过补集转化，问题转化为三个区间共同出现的数的数量。既然是膜队，就要对每个区间分别处理。现在我们需要维护的信息就十分明确了：维护资瓷快速求交集的区间出现元素。这用数据结构并不好维护，所以考虑万能的$\texttt{bitset}$。现在问题来了：$\texttt{bitset}$只能对于每个元素，维护是否出现过，并不能维护出现多少次。~~此时情况开始变得辣手起来。。。~~

能不能通过某种手段，使得$\texttt{bitset}$中留够足够的空间使得每个相同的数能够连续存储在一段空间中？珂以！我们只需要在离散化时不去重就珂以了！在存储时，令$x+cnt[x]-1$这一位为1，就珂以了。（我在代码中是采用$\texttt{shadowice1984}$大佬的方法，令离散化之后的值为数列中小于等于该数的数的个数，就令$x-cnt[x]+1$这一位为1。显然这两种方法是基本相同的。那么，我们就做完啦！

~~才不是呢！lxl的题怎么能这么轻易做完？~~

我们发现我们的最终结果是保存在$1e5个​$$\texttt{bitset}​$里面的。算一波空间，会发现我们开不下。。。所以，我们的算法错了么？考虑充足的时限，我们把询问分成三次处理。这里还有一个避免分类讨论的trick，就是首先把询问储存下来，每次查询时将莫队的询问数组清零，并将当前处理的区间复制到莫队的数组里。这样能使代码简洁许多，并提高珂复用性。

### 我的代码：

```cpp
#include <cstdio>
#include <cmath>
#include <cstring>
#include <bitset>
#include <map>
#include <algorithm>

using std::sort;
using std::lower_bound;
using std::bitset;
using std::map;
using std::min;

const int maxn=1e5+100;

int a[maxn],blo[maxn],col[maxn],len[maxn];
bitset<maxn> ans[40000],tmp;
int l1[maxn],l2[maxn],l3[maxn],r1[maxn],r2[maxn],r3[maxn];

struct Q
{
    int l,r,id;
    bool operator< (const Q& q) const
    {
        return blo[l]==blo[q.l]?r<q.r:l<q.l;
    }
}qry[maxn*3];

inline void upd(int p,int x)
{
    if (x==1)
        tmp[a[p]-(col[a[p]])]=1,++col[a[p]];
    else
        tmp[a[p]-col[a[p]]+1]=0,--col[a[p]];
}

inline void init(int l,int r,int &tmp)
{
    tmp=0;
    for (int i=l;i<=r;++i)
    {
        ++tmp;
        qry[tmp].l=l1[i];qry[tmp].r=r1[i];qry[tmp].id=i-l+1;
        ++tmp;
        qry[tmp].l=l2[i];qry[tmp].r=r2[i];qry[tmp].id=i-l+1;
        ++tmp;
        qry[tmp].l=l3[i];qry[tmp].r=r3[i];qry[tmp].id=i-l+1;
        len[i-l+1]=r1[i]-l1[i]+1+r2[i]-l2[i]+1+r3[i]-l3[i]+1;
    }
}

void solve(const int l,const int r)
{
    int tmp=0;
    int L=1,R=0;
    // int leftid=(left-1)*3+1,rightid=min(tmp,(right-1)*3+1);
    // sort(qry+left,qry+right+1);
    init(l,r,tmp);
    sort(qry+1,qry+tmp+1);
    memset(col,0,sizeof(col));
    for (int i=1;i<=r-l+1;++i) ans[i].set();
    ::tmp.reset();
    for (int i=1;i<=tmp;++i)
    {
        const int l=qry[i].l,r=qry[i].r;
        while (R<r) upd(++R,1);
        while (L>l) upd(--L,1);
        while (R>r) upd(R--,-1);
        while (L<l) upd(L++,-1);
        // if (cnt==1) ans[qry[i].id]&=tmp;
        // else if (cnt==2) ans[qry[i].id-33333]&=tmp;
        // else ans[qry[i].id-(33333+33336-1)]&=tmp;
        ans[qry[i].id]&=::tmp;
    }
    // for (int i=1;i<=tmp;++i)
    //     len[qry[i].id]+=qry[i].r-qry[i].l+1;
    for (int i=1;i<=r-l+1;++i)
        printf("%d\n",len[i]-(int)ans[i].count()*3);
}

int main()
{
    int n,m;
    scanf("%d%d",&n,&m);
    map<int,int> mp;
    int T=sqrt(n);
    for (int i=1;i<=n;++i)
        blo[i]=(i-1)/T+1;
    for (int i=1;i<=n;++i)
        scanf("%d",a+i),++mp[a[i]];
    for (map<int,int>::iterator it=mp.begin(),lst=it++;it!=mp.end();lst=it++)
        it->second+=lst->second;
    for (int i=1;i<=n;++i) a[i]=mp[a[i]];
    for (int i=1;i<=m;++i)
        scanf("%d%d%d%d%d%d",l1+i,r1+i,l2+i,r2+i,l3+i,r3+i);
    for (int i=1;i<=m;i+=36000) solve(i,min(i+36000-1,m));
}
```

