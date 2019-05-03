---
title: 蒟蒻_WA自动机的模板库
date: 2022-01-21 21:33:28
tags:
categories: 
    - 模板
comments:
---

# 更新日志

2019.2.16

* 修复珂朵莉树代码中的错误
* FFT板子更新为预处理单位复根的版本(多项式基本操作请移步多项式算法总结qwq)

2019.2.17

* 新增NTT板子

<!-- more -->
2019.2.25

* 新增替罪羊树板子

2019.3.6

* 新增K-D Tree(2-D Tree) -> [简单题AC代码]
* 更新高消板子
* 更新LCT板子

2019.3.31
* 新增Miller-Rabin素数判断
* 新增Pollard-Rho大数分解

2019.4.4
* 新增SAM板子

2019.4.22
* 新增二维凸包
* 新增笛卡尔树
* 新增Manacher

2019.4.26
* 新增毒瘤圆方树

2019.4.27
* 新增广义圆方树

# 数学

## 线性筛

```cpp
inline void sieve(int n)
{
	for (int i=2;i<=n;++i)
	{
		if (!v[i]) {v[i]=i;prime[++cnt]=i;}
		for (int j=1;j<=cnt;++j)
		{
			if (prime[j]>v[i] || i*prime[j]>n) break;
			v[i*prime[j]]=prime[j];
		}
	}
	for (int i=1;i<=cnt;++i)
		isprime[prime[i]]=true;
}
```


## 高斯消元 
> 模板题【SDOI2006】异或方程组
```cpp
#include <cstdio>
#include <cmath>
#include <algorithm>

using std::fabs;
using std::swap;

const int maxn=1e3+10;
const double eps=1e-6;

inline int Gauss_Elimination(double (*A)[maxn],double* f,int n)
{
    for (int i=1,c=1,j;i<=n;++i)
    {
        for (j=c;j<=n && fabs(A[j][i])<eps;++j);
        if (j==n+1) continue;
        for (int k=1;k<=n+1;++k) swap(A[c][k],A[j][k]);
        for (int j=c+1;j<=n;++j)
            if (fabs(A[j][i])>eps) 
            {
                double t=A[j][i]/A[c][i];
                for (int k=i;k<=n+1;++k)
                    A[j][k]-=t*A[c][k];
            }
        ++c;
    }
    bool NoAnswer=false,InfAnswer=false;
    for (int i=n;i;--i)
    {
        bool NoVariables=true;
        for (int j=i;j<=n;++j)
            if (fabs(A[i][j])>eps) NoVariables=false;
        if (NoVariables)
            if (fabs(A[i][n+1])>eps) NoAnswer=true; // 0=C,C!=0,无解
            else InfAnswer=true; // 0=0,无穷多组解
        else
        {
            for (int j=i+1;j<=n;++j) A[i][n+1]-=A[i][j]*f[j];
            f[i]=A[i][n+1]/A[i][i];
        }
    }
    if (NoAnswer) return -1; // 无解返回-1.. 
    return !InfAnswer; //无穷多解返回0，有唯一解返回1.
}

int main()
{
    static double A[maxn][maxn],f[maxn];
    int n;
    scanf("%d",&n);
    for (int i=1;i<=n;++i)
        for (int j=1;j<=n+1;++j)
            scanf("%lf",&A[i][j]);
    int result=Gauss_Elimination(A,f,n);
    if (result^1) return printf("%d\n",result)&0;
    for (int i=1;i<=n;++i) printf("x%d=%.2lf\n",i,f[i]);
}
```



## 三分

```cpp
inline double F(double x)
{
	double f=0;
	for (int i=n;~i;--i)
		f=f*x+a[i];
	return f;
}

int main()
{
	double l,r;
	scanf("%d%lf%lf",&n,&l,&r);
	for (int i=n;~i;--i)
		scanf("%lf",a+i);
	while (l+eps<r)
	{
		double m1=l+(r-l)/3;
		double m2=r-(r-l)/3;
		if (F(m1)>F(m2)) r=m2;
			else l=m1;
	}
	printf("%.5lf",l);
}
```



## 矩阵快速幂

```cpp
Matrix operator^ (ll k)
{
	Matrix ans(n,m);
	for (int i=1;i<=n;++i)
		ans.a[i][i]=1;
	Matrix t=*this;
	for (;k;k>>=1)
	{
		if (k&1) ans=ans*t;
		t=t*t;
	}
	return ans;
}
```



## 乘法逆元

### 线性递推

```cpp
inv[1]=1;
inv[i]=(p-p/i)*inv[p%i]%p;
```

###  阶乘逆元

$\text{inv}(i)=\text{inv}(i+1) \times(i+1)$

## 有理数取模

```cpp
inline ll pow(int a,int b,int mod)
{
	ll ans=1ll;
	for (;b;b>>=1)
	{
		if (b&1) ans=ans*a%mod;
		a=(ll)a*a%mod;
	}
	return ans;
}

inline ll read()
{
	char ch;
	while (!isdigit(ch=getchar()));
	ll x=ch-48;
	while (isdigit(ch=getchar()))
		x=(x*10+ch-48)%mod;
	return x;
}

int main()
{
	ll a,b;
	a=read();b=read();
	if (!b) return puts("Angry!"),0; 
	printf("%lld",(ll)a*pow(b,mod-2,mod)%mod);
}
```

## Miller-Rabin
```cpp
int pr[]={2,3,5,7,11,13,17,19,23,29,31,37};

inline ll qpow(ll a,ll b,ll p)
{
    ll ans=1%p;
    for (;b;b>>=1)
    {
        if (b&1) ans=ans*a%p;
        a=a*a%p;
    }
    return ans;
}

inline bool miller_rabin(int n)
{
    if (n==1) return false;
    for (int i=0;i<12;++i) if (n==pr[i]) return true;
    int m=(n-1),k=0;
    while (!(m&1)) m>>=1,++k;
    for (int i=0;i<12 && pr[i]<n;++i)
    {
        ll x=qpow(pr[i],m,n),y=x;
        for (int t=0;t<k;++t)
        {
            x=x*x%n;
            if (x==1 && y!=1 && y!=n-1) return false;
            y=x;
        }
        if (x!=1) return false;
    }
    return true;
}

```

## Pollard-Rho
```cpp
#include <cstdio>
#include <ctime>
#include <cmath>
#include <random>
#include <chrono>

using std::abs;

std::mt19937_64 Rnd(std::chrono::steady_clock::now().time_since_epoch().count());

typedef long long ll;
typedef unsigned long long ull;

int pr[]={2,3,5,7,11,13,17,19,23,29,31,37};

inline ll gcd(ll a,ll b){return b==0?a:gcd(b,a%b);}

inline ll rnd(ll x){return (ll)(Rnd()%x+1);}

inline ll quick_pow(ll a,ll b,ll p)
{
    ll ans=1%p;
    for (;b;b>>=1)
    {
        if (b&1) ans=ans*a%p;
        a=a*a%p;
    }
    return ans;
}

inline ll slow_mul(ull a,ll b,ll p)
{
    ull ans=0;
    for (;b;b>>=1)
    {
        if (b&1) ans=(ull)(ans+a);
        if (ans>=p) ans-=p;
        a=(ull)(a+a)%p;
        if (a>=p) a-=p;
    }
    return ans;
}

inline bool miller_rabin(ll n)
{
    if (n==1) return false;
    for (int i=0;i<12;++i) if (n==pr[i]) return true;
    if (n%2==0 || n%3==0 || n%5==0) return false;
    ll m=(n-1),k=0;
    while (!(m&1)) m>>=1,++k;
    for (int i=0;i<12 && pr[i]<n;++i)
    {
        ll x=quick_pow(pr[i],m,n),y=x;
        for (int t=0;t<k;++t)
        {
            x=slow_mul(x,x,n);
            if (x==1 && y!=1 && y!=n-1) return false;
            y=x;
        }
        if (x!=1) return false;
    }
    return true;
}

#define f(x) ((slow_mul(x,x,n)+c)%n)
inline ll rho(ll n)
{
    if (!(n&1)) return 2;
    if (n%3==0) return 3;
    ll x=0,y=0,t=1,c=rnd(n-1),q=1;
    for (int k=2;;k<<=1,y=x,q=1)
    {
        for (int i=1;i<k;++i)
        {
            x=f(x);
            q=slow_mul(q,abs(x-y),n);
            if (!(i&0x7f))
                if ((t=gcd(q,n))>1) break;
        }
        if (t>1 || (t=gcd(q,n))>1) break; 
    }
    return t;
}

ll max_p;

ll solve(ll n)
{
    if (n==1) return 1;
    if (miller_rabin(n)) return max_p=n>max_p?n:max_p; 
    ll t=n;
    while (t==n) t=rho(n);
    solve(n/t);solve(t);
    return max_p;
}

int main()
{
    freopen("pol.in","r",stdin);
    freopen("pol.out","w",stdout);
    int T;
    scanf("%d",&T);
    while (T--)
    {
        ll x;
        scanf("%lld",&x);max_p=0;
        ll ret=solve(x);
        if (ret==x) puts("Prime");
        else printf("%lld\n",ret);
    }
    fprintf(stderr,"%d",clock());
}
```

## FFT

```cpp
#include <cstdio>
#include <cmath>

const double Pi=acos(-1.0);
const int maxn=2e6+100;

double q[maxn];
int limit=1,rev[maxn];

struct Complex
{
    double real,imag;
    Complex(double real,double imag):real(real),imag(imag){}
    Complex(){}
    Complex conj();
}w[maxn],winv[maxn],A[maxn];

inline Complex Complex::conj(){return Complex(real,-imag);}
inline Complex operator+(const Complex& a,const Complex& b){return Complex(a.real+b.real,a.imag+b.imag);}
inline Complex operator-(const Complex& a,const Complex& b){return Complex(a.real-b.real,a.imag-b.imag);}
inline Complex operator*(const Complex& a,const Complex& b){return Complex(a.real*b.real-a.imag*b.imag,a.real*b.imag+a.imag*b.real);}

template<typename T>
inline void swap(T& a,T& b){T t=a;a=b;b=t;}

inline void DFT(Complex* A,Complex* w,int limit)
{
    for (int i=0;i<limit;++i)
        if (i<rev[i]) swap(A[i],A[rev[i]]);
    for (int mid=1;mid<limit;mid<<=1)
        for (int R=mid<<1,j=0;j<limit;j+=R)
            for (int k=0;k<mid;++k)
            {
                Complex x=A[j+k],y=w[limit/mid/2*k]*A[j+mid+k];
                A[j+k]=x+y;
                A[j+mid+k]=x-y;
            }
}

inline void prework(int n)
{
    int l=0;
    while (limit<=(n<<1)+1) limit<<=1,++l;
    for (int i=0;i<limit;++i)
        rev[i]=(rev[i>>1]>>1)|((i&1)<<(l-1));
    for (int i=0;i<limit;++i) 
        w[i]=Complex(cos(Pi*2/limit*i),sin(Pi*2/limit*i)),winv[i]=w[i].conj();
}

int main()
{
    int n,m;
    scanf("%d%d",&n,&m);
    for (int i=0;i<=n;++i)
        scanf("%lf",&A[i].real);
    for (int i=0;i<=m;++i)
        scanf("%lf",&A[i].imag);
    prework(n>=m?n:m);
    DFT(A,w,limit);
    for (int i=0;i<limit;++i)
        A[i]=A[i]*A[i];
    DFT(A,winv,limit);
    for (int i=0;i<=n+m;++i)
        printf("%d ",(int)(A[i].imag/2/limit+0.1));
}
```



## NTT

```cpp
inline int qpow(int a,int b,int p)
{
    int ans=1%p;
    for (;b;b>>=1,a=(ll)a*a%p)
        if (b&1) ans=(ll)ans*a%p;
    return ans;
}

inline void prework(int n,int m)
{
    int l=0;
    while (limit<=(n+m+1)) limit<<=1,++l;
    w[0]=1;w[1]=qpow(g,(P-1)/limit,P),winv[0]=1,winv[1]=qpow(w[1],P-2,P);
    for (int i=2;i<limit;++i)
        w[i]=(ll)w[i-1]*w[1]%P,winv[i]=1ll*winv[i-1]*winv[1]%P;
    for (int i=1;i<limit;++i)
        rev[i]=(rev[i>>1]>>1)|((i&1)<<(l-1));
}

inline void NTT(int *A,int *w,int limit)
{
    for (int i=0;i<limit;++i)
        if (i<rev[i]) swap(A[i],A[rev[i]]);
    for (int mid=1;mid<limit;mid<<=1)
        for (int R=mid<<1,j=0;j<limit;j+=R)
            for (int k=0;k<mid;++k)
            {
                int x=A[j+k],y=(ll)A[j+k+mid]*w[limit/2/mid*k]%P;
                A[j+k]=(x+y)%P;A[j+mid+k]=(x-y+P)%P;
            }
}

inline void DFT(int *A){ NTT(A,w,limit); }

inline void IDFT(int *A)
{
    NTT(A,winv,limit);
    int inv=qpow(limit,P-2,P);
    for (int i=0;i<=limit;++i)
        A[i]=((ll)A[i]*inv)%P;
}
```

# 计算几何
## 二维凸包
```cpp
#include <cstdio>
#include <cmath>
#include <algorithm>

using std::sort;
using std::sqrt;

const int maxn=1e5+100;

struct Point
{
    double x,y;
    Point(double x,double y):x(x),y(y){}
    Point(){}
    bool operator< (const Point& p) const{return x==p.x?y<p.y:x<p.x;}
    Point operator- (const Point& p){return Point(x-p.x,y-p.y);}
    double operator* (const Point& p){return x*p.x+y*p.y;}
}p[maxn];

typedef Point Vector;

inline double dis(const Point& a,const Point& b){return sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y));}

inline double cross(const Vector& a,const Vector& b){return a.x*b.y-b.x*a.y;}

int main()
{
    int n;
    scanf("%d",&n);
    for (int i=1;i<=n;++i)
        scanf("%lf%lf",&p[i].x,&p[i].y);
    sort(p+1,p+n+1);
    static int stack[maxn];
    int top=0;
    stack[++top]=1;
    bool used[maxn];
    for (int i=2;i<=n;++i)
    {
        while (top>1 && cross(p[stack[top]]-p[stack[top-1]],p[i]-p[stack[top]])<=0) used[stack[top--]]=false;
        stack[++top]=i;
        used[i]=true;
    }
    int prec=top;
    for (int i=n;i;--i)
        if (!used[i])
        {
            while (top>prec && cross(p[stack[top]]-p[stack[top-1]],p[i]-p[stack[top]])<=0) used[stack[top--]]=false;
            stack[++top]=i;
            used[i]=true;
        }
    double ans=0;
    for (int i=2;i<=top;++i)
        ans+=dis(p[stack[i]],p[stack[i-1]]);
    printf("%.2lf",ans);
}
```

# 字符串

## Manacher
```cpp
#include <cstdio>
#include <cstring>

const int maxn=3e7+10;

char buf[maxn],s[maxn];
int f[maxn];

template<class T>inline T min(T a,T b){return a<b?a:b;}
template<class T>inline T max(T a,T b){return a<b?b:a;}

inline void manacher(int n)
{
    int maxright=0,mid=0;
    for (int i=1;i<n;++i)
    {
        if (i<maxright) f[i]=min(f[(mid<<1)-i],f[mid]+mid-i);
        else f[i]=1;
        while (s[i+f[i]]==s[i-f[i]]) ++f[i];
        if (i+f[i]>maxright) maxright=i+f[i],mid=i;
    }
}

int main()
{
    int n;
    scanf("%s",buf);
    n=strlen(buf);
    s[0]=s[1]='@';
    for (int i=0;i<n;++i)
        s[(i<<1)+2]=buf[i],s[(i<<1)+3]='@';
    n=n*2+2;s[n]=0;
    manacher(n);
    int ans=0;
    for (int i=0;i<n;++i)
        ans=max(ans,f[i]);
    printf("%d",ans-1);
}
```

## 制胡窜哈希

```cpp
inline void hs1(char* s)
{
	pw1[0]=1;hsh1[0]=s[0];
	for (int i=1;i<n;++i)
	{
		pw1[i]=(ll)pw1[i-1]*seed%mod1;
		hsh1[i]=((ll)hsh1[i-1]*seed+s[i])%mod1;
	}
}

inline void hs2(char* s)
{
	pw2[0]=1;hsh2[0]=s[0];
	for (int i=1;i<n;++i)
	{
		pw2[i]=(ll)pw2[i-1]*seed%mod2;
		hsh2[i]=((ll)hsh2[i-1]*seed+s[i])%mod2;
	}
}

inline pair<int,int> gethash(int l,int r)
{
	int t1=((hsh1[r]-(ll)hsh1[l-1]*pw1[r-l+1]%mod1)+mod1)%mod1;
	int t2=((hsh2[r]-(ll)hsh2[l-1]*pw2[r-l+1]%mod2)+mod2)%mod2;
	return make_pair(t1,t2);
}
```



## KMP

```cpp
inline void getfail(int n)
{
    for (int i=1,j;i<n;++i)
    {
        j=f[i];
        while (j && P[j]!=P[i]) j=f[j];
        f[i+1]=P[i]==P[j]?j+1:0;
    }
}

inline void kmp(int n,int m)
{
    for (int i=0,j=0;i<n;++i)
    {
        while (j && T[i]!=P[j]) j=f[j];
        if (T[i]==P[j]) ++j;
        if (j==m) printf("%d\n",i-m+2);
    }
}
```



## AC自动机

```cpp
void add(const char* s)
{
	int n=strlen(s),u=0;
	for (int i=0;i<n;++i)
	{
		int c=idx(s[i]);
		if (!ch[u][c]) ch[u][c]=++cnt;
		u=ch[u][c];
	}
	++tag[u];
}

void getfail()
{
	queue<int> q;
	for (int i=0;i<26;++i)
		if (ch[0][i]) q.push(ch[0][i]);
	while (!q.empty())
	{
		int u=q.front();q.pop();
		for (int i=0;i<26;++i)
		{
			int c=ch[u][i];
			if (!c) {ch[u][i]=ch[f[u]][i];continue;}
			q.push(c);
			int v=f[u];
			while (v && !ch[v][i]) v=f[v];
			f[c]=ch[v][i];
			last[c]=tag[f[c]]?f[c]:last[f[c]];
		}
	}
}

int get(const char *s)
{
    int n=strlen(s);
    int x=0,ret=0;
    for(int i=0;i<n;++i)
    {
        x=ch[x][idx(s[i])];
        for(int j=x;j;j=last[j])
        	if(!vis[j])vis[j]=1,ret+=tag[j];
    }
    return ret;
}
```



## SA

```cpp
inline void build_sa(int n,int m)
{
    int *x=t1,*y=t2;
    for (int i=0;i<m;++i) c[i]=0;
    for (int i=0;i<n;++i) c[x[i]=s[i]]++;
    for (int i=1;i<m;++i) c[i]+=c[i-1];
    for (int i=n-1;~i;--i) sa[--c[x[i]]]=i;
    for (int k=1,p=1;k<=n && p<n;k<<=1,m=p)
    {
        p=0;
        for (int i=n-k;i<n;++i) y[p++]=i;
        for (int i=0;i<n;++i) if (sa[i]>=k) y[p++]=sa[i]-k;
        for (int i=0;i<m;++i) c[i]=0;
        for (int i=0;i<n;++i) c[t3[i]=x[y[i]]]++;
        for (int i=1;i<m;++i) c[i]+=c[i-1];
        for (int i=n-1;~i;--i) sa[--c[t3[i]]]=y[i];
        p=1;swap(x,y);x[sa[0]]=0;
        for (register int i=1;i<n;++i)
            x[sa[i]]=(y[sa[i]+k]==y[sa[i-1]+k] && y[sa[i]]==y[sa[i-1]])?p-1:p++;
    }
}

inline void get_height(int n)
{
    int k=0,j=0;
    for (int i=0;i<=n;++i) rank[sa[i]]=i;
    for (int i=0;i<n;height[rank[i++]]=k)
        for (k?--k:0,j=sa[rank[i]-1];s[j+k]==s[i+k];++k);
}
```

## SAM
```cpp
#include <cstdio>
#include <cstring>

const int maxn=2e6+1000;

int tr[maxn][26],parent[maxn],mx[maxn],right[maxn],cnt=1,last=1;

inline void radixsort(int n)
{
    static int c[maxn],id[maxn];
    for (int i=1;i<=cnt;++i) ++c[mx[i]];
    for (int i=1;i<=n;++i) c[i]+=c[i-1];
    for (int i=cnt;i;--i) id[--c[mx[i]]]=i;
    for (int i=cnt-1;~i;--i) right[parent[id[i]]]+=right[id[i]];
}

inline void insert(int x)
{
    int p=last,np=last=++cnt;
    right[np]=1;mx[np]=mx[p]+1;
    while (p && !tr[p][x]) tr[p][x]=np,p=parent[p];
    if (!p) parent[np]=1;
    else
    {
        int q=tr[p][x];
        if (mx[q]==mx[p]+1) parent[np]=q;
        else
        {
            int nq=++cnt;
            mx[nq]=mx[p]+1;
            memcpy(tr[nq],tr[q],sizeof(tr[q]));
            while (p && tr[p][x]==q) tr[p][x]=nq,p=parent[p];
            parent[nq]=parent[q];parent[q]=parent[np]=nq;
        }
    }
}

int main()
{
    int n;
    static char s[maxn];
    scanf("%s",s);
    n=strlen(s);
    for (int i=0;i<n;++i) insert(s[i]-'a');
    radixsort(n);
    int ans=0;
    for (int i=1;i<=cnt;++i) if (right[i]>1 && right[i]*mx[i]>ans) ans=right[i]*mx[i]; 
    printf("%d\n",ans);
}
```

# 图论
## 广义圆方树(APIO2018 铁人两项)
```cpp
#include <cstdio>
#include <vector>
#include <algorithm>

using std::vector;

const int maxn=2e6+100;

vector<int> G[maxn],T[maxn];
int dfn[maxn],low[maxn],dfc,tot,cnt,val[maxn],siz[maxn];

template<class T>inline T max(T a,T b){return a<b?b:a;}
template<class T>inline T min(T a,T b){return a<b?a:b;}
template<class T>inline void swap(T& a,T& b){a^=b^=a^=b;}

inline void tarjan(int u,int ff)
{
    static int stack[maxn],top=0;
    dfn[u]=low[u]=++dfc;
    stack[++top]=u;
    ++cnt;val[u]=-1;
    for (auto v:G[u])
    {
        if (v==ff) continue;
        if (!dfn[v])
        {
            tarjan(v,u);
            low[u]=min(low[u],low[v]);
            if (low[v]>=dfn[u])
            {
                ++tot;
                int x;
                do
                {
                    x=stack[top--];
                    ++val[tot];
                    T[x].push_back(tot);
                    T[tot].push_back(x);
                } while (x!=v);
                T[u].push_back(tot);
                T[tot].push_back(u);
                ++val[tot];
            }
        }
        else low[u]=min(low[u],dfn[v]);
    }
}

inline int64_t calc(int u,int ff,int n,int N)
{
    int64_t ans=0;
    siz[u]=(u<=N);
    for (auto v:T[u])
    {
        if (v==ff) continue;
        ans+=calc(v,u,n,N);
        ans+=(int64_t)siz[u]*siz[v]*val[u];
        siz[u]+=siz[v];
    }
    ans+=(int64_t)siz[u]*(n-siz[u])*val[u];
    return ans;
}

int main()
{
    int n,m;
    scanf("%d%d",&n,&m);
    tot=n;
    for (int i=1,u,v;i<=m;++i)
        scanf("%d%d",&u,&v),G[u].push_back(v),G[v].push_back(u);
    int64_t ans=0;
    for (int i=1;i<=n;++i)
        if (!dfn[i])
        {
            cnt=0;
            tarjan(i,0);
            ans+=calc(i,0,cnt,n);
        } 
    printf("%lld\n",ans*2);
}
```

## 静态仙人掌（圆方树）
```cpp
#include <cstdio>
#include <algorithm>

const int maxn=3e4+1200;

struct Graph
{
    struct Edge
    {
        int to,next,w;
    }edge[maxn<<1];
    int head[maxn],cnt;

    inline void _add(int u,int v,int w)
    {
        edge[++cnt].next=head[u];
        edge[cnt].to=v;
        edge[cnt].w=w;
        head[u]=cnt;
    }
    inline void add(int u,int v,int w){_add(u,v,w);_add(v,u,w);}
}G,T;

int dfn[maxn],low[maxn],tot,dfc;
int fa[maxn],val[maxn],sum[maxn];

template<class T>inline T min(T a,T b){return a<b?a:b;}
template<class T>inline T max(T a,T b){return a<b?b:a;}
template<class T>inline void swap(T &a,T &b){a^=b^=a^=b;}

inline void build(int u,int v,int w)
{
    ++tot;
    int s=w;
    for (int x=v;x!=fa[u];x=fa[x])
        sum[x]=s,s+=val[x];
    sum[tot]=sum[u];
    sum[u]=0;
    for (int x=v;x!=fa[u];x=fa[x])
    {
        int sp=min(sum[x],sum[tot]-sum[x]);
        T.add(x,tot,sp);
    }
}

inline void tarjan(int u,int ff)
{
    dfn[u]=low[u]=++dfc;
    for (int i=G.head[u];i;i=G.edge[i].next)
    {
        int v=G.edge[i].to;
        if (v==ff) continue;
        if (!dfn[v])
        {
            fa[v]=u;
            val[v]=G.edge[i].w;
            tarjan(v,u);
            low[u]=min(low[u],low[v]);
        }
        else low[u]=min(low[u],dfn[v]);
        if (low[v]<=dfn[u]) continue;
        T.add(u,v,G.edge[i].w);
    }
    for (int i=G.head[u];i;i=G.edge[i].next)
    {
        int v=G.edge[i].to;
        if (fa[v]!=u && dfn[v]>dfn[u])
            build(u,v,G.edge[i].w);
    }
}

int depth[maxn],siz[maxn],top[maxn],son[maxn],dis[maxn];

inline void dfs(int u,int ff,int dep,int di)
{
    depth[u]=dep;
    dis[u]=di;
    fa[u]=ff;
    int maxs=-1;
    siz[u]=1;
    for (int i=T.head[u];i;i=T.edge[i].next)
    {
        int v=T.edge[i].to;
        if (v!=ff)
        {
            dfs(v,u,dep+1,di+T.edge[i].w);
            siz[u]+=siz[v];
            if (siz[v]>maxs) maxs=siz[v],son[u]=v;
        }
    }
}

inline void dfs(int u,int topf)
{
    top[u]=topf;
    if (!son[u]) return;
    dfs(son[u],topf);
    for (int i=T.head[u];i;i=T.edge[i].next)
    {
        int v=T.edge[i].to;
        if (v!=fa[u] && v!=son[u])
            dfs(v,v);
    }
}

inline int LCA(int u,int v)
{
    while (top[u]!=top[v])
    {
        if (depth[top[u]]<depth[top[v]]) swap(u,v);
        u=fa[top[u]];
    }
    return depth[u]>depth[v]?v:u;
}

inline int find(int u,int lca)
{
    int res=son[lca];
    while (top[u]!=top[lca])
        res=top[u],u=fa[top[u]];
    return u==lca?res:son[lca];
}

inline int query(int u,int v,int n)
{
    int lca=LCA(u,v);
    if (lca<=n) return dis[u]+dis[v]-2*dis[lca];
    int x=find(u,lca),y=find(v,lca);
    int ans=dis[u]+dis[v]-dis[x]-dis[y];
    if (sum[x]<sum[y]) swap(x,y);
    ans+=min(sum[x]-sum[y],sum[lca]-sum[x]+sum[y]);
    return ans;
}

int main()
{
    int n,m,q;
    scanf("%d%d%d",&n,&m,&q);
    for (int i=1,u,v,w;i<=m;++i)
        scanf("%d%d%d",&u,&v,&w),G.add(u,v,w);
    tot=n;
    tarjan(1,0);
    dfs(1,0,1,0);
    dfs(1,1);
    for (int i=1,u,v;i<=q;++i)
    {
        scanf("%d%d",&u,&v);
        printf("%d\n",query(u,v,n));
    }
}
```

## Kruskal

```cpp
struct Edge
{
	int u,v,w;
	bool operator< (const Edge& e) const {return w<e.w;}
}edge[maxn];
int fa[maxn],cnt;

inline int find(int x){return x==fa[x]?x:fa[x]=find(fa[x]);}

inline int kruskal(int m)
{
	int ans=0;
	std::sort(edge+1,edge+m+1);
	for (int i=1;i<=m;++i)
	{
		int u=edge[i].u,v=edge[i].v,w=edge[i].w;
		if (find(u)!=find(v))
			++cnt,fa[find(u)]=find(v),ans+=w;
	}
	return ans;
}
```

## LCA

###  树剖  

```cpp
void dfs(int u,int f,int dep)
{
    depth[u]=dep;
    fa[u]=f;
    siz[u]=1;
    int maxs=-1;
    for (int i=head[u];i;i=edge[i].next)
    {
        int v=edge[i].to;
        if (v!=f) 
        {
            dfs(v,u,dep+1);
            siz[u]+=siz[v];
            if (siz[v]>=maxs) maxs=siz[v],son[u]=v;
        }
    }
}

void dfs(int u,int topf)
{
    top[u]=topf;
    if (!son[u]) return;
    dfs(son[u],topf);
    for (int i=head[u];i;i=edge[i].next)
    {
        int v=edge[i].to;
        if (v!=fa[u] && v!=son[u]) dfs(v,v);
    }
}

int lca(int u,int v)
{
    while (top[u]!=top[v])
    {
        if (depth[top[u]]<depth[top[v]]) swap(u,v);
        u=fa[top[u]];
    }
    if (depth[u]>depth[v]) swap(u,v);
    return u;
}
```

### 倍增

```cpp
inline void dfs(int u,int fa,int dep)
{
    f[u][0]=fa;
    depth[u]=dep;
    for (int i=1;(1<<i)<=n;++i)
        f[u][i]=f[f[u][i-1]][i-1];
    for (int i=head[u];i;i=edge[i].next)
    {
        int v=edge[i].to;
        if (v!=fa) dfs(v,u,dep+1);
    }
}

inline int lca(int u,int v)
{
    if (depth[u]<depth[v]) swap(u,v);
    int l=depth[u]-depth[v];
    for (int i=0;(1<<i)<=l;++i)
        if (l&(1<<i)) u=f[u][i];
    if (u==v) return u;
    for (int i=22;i>=0;--i)
        if (f[u][i]!=f[v][i]) 
            u=f[u][i],v=f[v][i];
    return f[u][0];
}
```

### DFS序转RMQ

```cpp
inline int Min(int x,int y){return depth[x]<depth[y]?x:y;}

void dfs(int u,int f,int dep)
{
    fa[u]=f;
    dfn[++tot]=u;
    pre[u]=tot;
    depth[tot]=dep;
    for (int i=head[u];i;i=edge[i].next)
    {
        int v=edge[i].to;
        if (v!=f)
        {
            dfs(v,u,dep+1);
            dfn[++tot]=u;
            depth[tot]=dep;
        }
    }
}

inline void rmq_init()
{
    for (int i=1;i<=tot;++i) st[i][0]=i; 
    for (int j=1;(1<<j)<=tot;++j)
        for (int i=1;i+(1<<j)-1<=tot;++i)
            st[i][j]=Min(st[i][j-1],st[i+(1<<(j-1))][j-1]);
}

inline int rmq(int L,int R)
{
    int k=0;
    while ((1<<(k+1))<=R-L+1) ++k;
    return Min(st[L][k],st[R-(1<<k)+1][k]);
}

inline int lca(int x,int y)
{
    x=pre[x];y=pre[y];
    if (x>y) swap(x,y);
    return dfn[rmq(x,y)];
}
```

### Tarjan

```cpp
inline int find(int x){ return x==fa[x]?x:fa[x]=find(fa[x]);}

void dfs(int u,int f)
{
    fa[u]=u;
    for (int i=head[u];i;i=edge[i].next)
    {
        int v=edge[i].to;
        if (v!=f) dfs(v,u),fa[v]=u;
    }
    for (int i=query_head[u];~i;i=query[i].next)
    {
        int v=query[i].v;
        if (vis[v]) query[i].ans=query[i^1].ans=find(v);
    } 
    vis[u]=true;
}
```



## 单源最短路径

> Dijkstra

```cpp
void Dijkstra(int s)
{
	memset(d,0x3f,sizeof(d));
	std::priority_queue<Node> q;
	d[s]=0;
	q.push(Node(s,0));
	while (!q.empty())
	{
		int u=q.top().u,di=q.top().dis;q.pop();
		if (di!=d[u]) continue;
		for (int i=head[u];i;i=edge[i].next)
		{
			int v=edge[i].to;
			if (d[v]>d[u]+edge[i].w)
			{
				d[v]=d[u]+edge[i].w;
				q.push(Node(v,d[v]));
			}
		}
	}
}
```

> SPFA

```cpp
void spfa()
{
    d[s]=0;
    vis[s]=1;
    queue<int> q;
    q.push(s);
    while (!q.empty())
    {
        int u=q.front(); q.pop();vis[u]=false;
        for (int i=head[u];i;i=edges[i].next)
        {
            Edge& e=edges[i];
            if (d[e.to]>d[u]+e.w)
            {
                d[e.to]=d[u]+e.w;
                if (!vis[e.to]) q.push(e.to),vis[e.to]=true;
            }
        }
    }
}
```



## 网络流

> Dinic最大流

```cpp
bool bfs()
{
	queue<int> q;
	q.push(s);
	memset(depth,0,sizeof(depth));
	depth[s]=1;
	while (!q.empty())
	{
		int u=q.front();q.pop();
		for (int i=head[u];~i;i=edge[i].next)
		{
			int v=edge[i].to;
			if (!depth[v] && edge[i].w>0)
				depth[v]=depth[u]+1,q.push(v);
		}
	}
	return depth[t]>0; 
}

int dfs(int u,int flow)
{
	if (u==t) return flow;
	for (int &i=cur[u];~i;i=edge[i].next)
	{
		int v=edge[i].to;
		if (depth[v]==depth[u]+1 && edge[i].w>0)
		{
			int d=dfs(v,min(flow,edge[i].w));
			if (d>0)
			{
				edge[i].w-=d;
				edge[i^1].w+=d;
				return d;
			}
		}
	}
	return 0;
}

int dinic(int n)
{
	int ans=0;
	while (bfs())
	{
		memcpy(cur,head,sizeof(cur));
		while (int d=dfs(s,INF)) ans+=d;
	}
	return ans;
}
```

> Edmonds-Karp费用流

```cpp
inline bool spfa()
{
	static int vis[maxn];
	queue<int> q;
	q.push(s);
	memset(vis,0,sizeof(vis));
	memset(d,0x3f,sizeof(d));
	d[s]=0;vis[s]=true;
	while (!q.empty())
	{
		int u=q.front();q.pop();
		vis[u]=false;
		for (int i=head[u];~i;i=edge[i].next)
		{
			int v=edge[i].to;
			if (edge[i].flow>0 && d[v]>d[u]+edge[i].cost)
			{
				d[v]=d[u]+edge[i].cost;
				incf[v]=edge[i].flow;
				inedge[v]=i;
				pre[v]=u;
				if (!vis[v]) vis[v]=true,q.push(v);
			}
		}
	}
	return d[t]<inf;
}

inline void mcmf()
{
	int flow=0,cost=0;
	while (spfa())
	{
		int u=t,mn=inf;
		for (int i=t;i!=s;i=pre[i])
			mn=min(mn,incf[i]);
		for (int i=t;i!=s;i=pre[i])
		{
			edge[inedge[i]].flow-=mn;
			edge[inedge[i]^1].flow+=mn;
		}
		cost+=d[t]*mn;
		flow+=mn;
	}
	printf("%d %d",flow,cost);
}

```



## 割点

```cpp
inline void tarjan(int u,int rt)
{
	int chd=0;
	dfn[u]=low[u]=++idx;
	for (int i=head[u];i;i=edge[i].next)
	{
		int v=edge[i].to;
		if (!dfn[v])
		{
			tarjan(v,rt);
			low[u]=min(low[u],low[v]);
			if (low[v]>=dfn[u])
			{
				++chd;
				if (chd>1 || u!=rt)
					iscut[u]=true;
			}
		}
		else low[u]=min(low[u],dfn[v]);
	}
}
```



## Tarjan缩点

```cpp
struct Graph
{
	int head[maxn],cnt;
	Edge edge[maxm];
	inline void add(int u,int v)
	{
		edge[++cnt].next=head[u];
		edge[cnt].to=v;
		head[u]=cnt;
	}
}old,dag;

void tarjan(int u)
{
	dfn[u]=low[u]=++idx;
	stk[++top]=u;
	instack[u]=true;
	for (int i=old.head[u];i;i=old.edge[i].next)
	{
		int v=old.edge[i].to;
		if (!dfn[v])
		{
			tarjan(v);
			low[u]=min(low[u],low[v]);
		}
		else if (instack[v]) low[u]=min(low[u],dfn[v]);
	}
	if (dfn[u]==low[u])
	{
		int v;
		++scc_cnt;
		do
		{
			v=stk[top--];
			belong[v]=scc_cnt;
			w[scc_cnt]+=val[v];
			instack[v]=false;
		}while (v!=u);
	}
}

inline void rebuild(int n)
{
	for (int u=1;u<=n;++u)
		for (int i=old.head[u];i;i=old.edge[i].next)
			if (belong[u]!=belong[old.edge[i].to])
				dag.add(belong[u],belong[old.edge[i].to]);
}

inline void work(int n)
{
	for (int i=1;i<=n;++i)
		if (!dfn[i]) tarjan(i);
}
```



## 2-SAT

```cpp
inline void tarjan(int u)
{
	dfn[u]=low[u]=++idx;
	instack[u]=true;
	stk[++top]=u;
	for (int i=head[u];i;i=edge[i].next)
	{
		int v=edge[i].to;
		if (instack[v]) low[u]=min(low[u],dfn[v]);
			else if (!dfn[v])
				tarjan(v),low[u]=min(low[u],low[v]);
	}
	if (dfn[u]==low[u])
	{
		int v;++scc_cnt;
		do
		{
			v=stk[top--];
			belong[v]=scc_cnt;
			instack[v]=false;
		}while (v!=u);
	}
}

inline bool twosat(int n)
{
	for (int i=1;i<=n<<1;++i)
		if (!dfn[i]) tarjan(i);
	for (int i=1;i<=n;++i)
		if (belong[i]==belong[i+n]) return false;
	return true;
}
```



# 数据结构

## ST表

```cpp
inline void prework(int n)
{
	for (rint j=1;(1<<j)<=n;++j)
		for (rint i=1;(i+(1<<j)-1)<=n;++i)
			d[i][j]=max(d[i][j-1],d[i+(1<<(j-1))][j-1]);
}

inline int query(int i,int j)
{
	int k=0,len=j-i+1;
	while (1<<(k+1)<=len)++k;
	return max(d[i][k],d[j-(1<<k)+1][k]);
}
```

 

## 线段树2

```cpp
void build(int l,int r,int o)
{
	mul[o]=1;
	if (l==r)
	{
		scanf("%lld",sumv+o);
		return;
	}
	int m=(l+r)>>1;
	build(ls);build(rs);
	pushup(o);
}

inline void pushdown(int o,int len)
{
	sumv[o<<1]=(sumv[o<<1]*mul[o]+addv[o]*(len-(len>>1)))%p;
	sumv[o<<1|1]=(sumv[o<<1|1]*mul[o]+addv[o]*(len>>1))%p;
	mul[o<<1]=mul[o<<1]*mul[o]%p;
	mul[o<<1|1]=mul[o<<1|1]*mul[o]%p;
	addv[o<<1]=(addv[o<<1]*mul[o]+addv[o])%p;
	addv[o<<1|1]=(addv[o<<1|1]*mul[o]+addv[o])%p;
	addv[o]=0;mul[o]=1;
}

void multiplicate(int L,int R,int c,int l,int r,int o)
{
	if (L<=l && R>=r)
	{
		mul[o]=mul[o]*c%p;
		addv[o]=addv[o]*c%p;
		sumv[o]=sumv[o]*c%p;
		return;
	}
	pushdown(o,r-l+1);
	int m=(l+r)>>1;
	if (L<=m) multiplicate(L,R,c,ls);
	if (R>m)  multiplicate(L,R,c,rs);
	pushup(o);
}

void add(int L,int R,int c,int l,int r,int o)
{
	if (L<=l && R>=r)
	{
		addv[o]=(addv[o]+c)%p;
		sumv[o]=(sumv[o]+c*(r-l+1))%p;
		return;
	}
	pushdown(o,r-l+1);
	int m=(l+r)>>1;
	if (L<=m) add(L,R,c,ls);
	if (R>m)  add(L,R,c,rs);
	pushup(o);
}

long long Querysum(int L,int R,int l,int r,int o)
{
	if (L<=l && R>=r) return sumv[o]%p;
	pushdown(o,r-l+1);
	int m=(l+r)>>1;
	long long tot=0;
	if (L<=m) tot=(tot+Querysum(L,R,ls))%p;
	if (R>m)  tot=(tot+Querysum(L,R,rs))%p;
	pushup(o);
	return tot;
}
```



## 左偏树

```cpp
inline int merge(int x,int y)
{
	if (!x || !y) return x+y;
	if (val[x]>val[y] || (val[x]==val[y] && x>y))
		swap(x,y);
	rs(x)=merge(rs(x),y);
	fa[rs(x)]=x;
	if (dis[rs(x)]>dis[ls(x)]) swap(ls(x),rs(x));
	dis[x]=dis[rs(x)]+1;
	return x;
}

inline int findroot(int x)
{
	while (fa[x]) x=fa[x];
	return x;
}

inline int pop(int x)
{
	int ret=val[x];
	fa[ls(x)]=fa[rs(x)]=0;
	val[x]=-1;
	merge(ls(x),rs(x));
	ls(x)=rs(x)=0;
	return ret;
}
```



## 主席树

```cpp
inline void insert(int x,int &rt,int oldrt,int l,int r)
{
    rt=++num;
    tree[rt]=tree[oldrt];
    ++tree[rt].sumv;
    if (l==r) return;
    int m=(l+r)>>1;
    if (x<=m) insert(x,tree[rt].ls,tree[oldrt].ls,l,m);
        else insert(x,tree[rt].rs,tree[oldrt].rs,m+1,r);
}

inline int query(int x,int lrt,int rrt,int l,int r)
{
    if (l==r) return l;
    int m=(l+r)>>1;
    int k=tree[tree[rrt].ls].sumv-tree[tree[lrt].ls].sumv;
    if (x<=k) return query(x,tree[lrt].ls,tree[rrt].ls,l,m);
        else return query(x-k,tree[lrt].rs,tree[rrt].rs,m+1,r);
}

```



## CDQ分治（三维偏序）

```cpp
#include <cstdio>
#include <algorithm>
#include <cstring>

using std::sort;

const int maxn=1e5+1000;

struct Tuple
{
    int a,b,c,cnt,ans;
    bool operator< (const Tuple& Tp) const 
    {
        if (a!=Tp.a) return a<Tp.a;
        if (b!=Tp.b) return b<Tp.b;
        return c<Tp.c;
    }
    bool operator!= (const Tuple& Tp)
    {
        return a!=Tp.a || b!=Tp.b || c!=Tp.c;
    }
}tmp[maxn],a[maxn];

struct cmp
{
    bool operator() (const Tuple& a,const Tuple& b)
    {
        if (a.b!=b.b) return a.b<b.b;
        return a.c<b.c;
    }
};

int c[maxn*2],n,k;

inline void update(int x,int y)
{
    for (int i=x;i<=k;i+=i&-i) c[i]+=y;
}

inline int query(int x)
{
    int ans=0;
    for (int i=x;i;i-=i&-i) ans+=c[i];
    return ans;
}

void solve(int l,int r)
{
    if (l==r) return;
    int mid=(l+r)>>1;
    solve(l,mid);solve(mid+1,r);
    sort(a+l,a+mid+1,cmp());sort(a+mid+1,a+r+1,cmp());
    for (int t1=l,t2=mid+1;t2<=r;++t2)
    {
        while (a[t1].b<=a[t2].b && t1<=mid) update(a[t1].c,a[t1].cnt),++t1;
        a[t2].ans+=query(a[t2].c);
    }
    for (int i=l;i<=mid;++i) 
        if (a[i].b<=a[r].b) update(a[i].c,-a[i].cnt);
            else break;
    // memset(c,0,sizeof(c));
}

int main()
{
    scanf("%d%d",&n,&k);
    for (int i=1;i<=n;++i) 
        scanf("%d%d%d",&tmp[i].a,&tmp[i].b,&tmp[i].c);
    std::sort(tmp+1,tmp+n+1);
    int tot;a[tot=1]=tmp[1];a[1].cnt=1;
    for (int i=2;i<=n;++i)
        if (tmp[i]!=tmp[i-1]) a[++tot]=tmp[i],a[tot].cnt=1;
            else ++a[tot].cnt;
    solve(1,tot);
    // for (int i=1;i<=tot;++i) printf("%d ",a[i].ans);
    static int ans[maxn];
    for (int i=1;i<=tot;++i) ans[a[i].ans+a[i].cnt-1]+=a[i].cnt;
    for (int i=0;i<n;++i) printf("%d\n",ans[i]);
}
```



## 点分治

```cpp
void getroot(int u,int fa)
{
    mxsiz[u]=0;siz[u]=1;
    for (int i=head[u];i;i=edge[i].next)
    {
        int v=edge[i].to;
        if (vis[v] || v==fa) continue;
        getroot(v,u);
        siz[u]+=siz[v];
        mxsiz[u]=max(mxsiz[u],siz[v]);
    }
    mxsiz[u]=max(mxsiz[u],S-siz[u]);
    if (mxsiz[u]<mxsiz[root]) root=u;
}

void getdis(int u,int fa,int d)
{
    tmp[++cnt]=d;
    for (int i=head[u];i;i=edge[i].next)
    {
        int v=edge[i].to;
        if (vis[v] || v==fa) continue;
        getdis(v,u,d+edge[i].w);
    }
}

void getans(int opt)
{
    sort(tmp+1,tmp+cnt+1);
    for (int u=1;u<=cnt;++u)
    {
        for (int t=1;t<=m;++t)
        {
            pair<int*,int*> p=equal_range(tmp+1,tmp+cnt+1,query[t]-tmp[u]);
            if (p.second!=p.first) count[t]+=opt*(p.second-p.first);
            // assert(p.second==p.first);
        }
    }
}

void solve(int u)
{
    vis[u]=1;
    getdis(u,cnt=0,0);
    getans(1);
    for (int i=head[u];i;i=edge[i].next)
    {
        int v=edge[i].to;
        if (vis[v]) continue;
        getdis(v,cnt=0,edge[i].w);
        getans(-1);
        mxsiz[root=0]=0x3f3f3f3f;
        S=siz[v];
        getroot(v,0);
        solve(v);
    }
}
```

## 笛卡尔树(Luogu P3793)
```cpp
#include <cstdio>

typedef long long unsigned int uint64_t;

namespace GenHelper
{
    unsigned z1,z2,z3,z4,b;
    unsigned rand_()
    {
        b=((z1<<6)^z1)>>13;
        z1=((z1&4294967294U)<<18)^b;
        b=((z2<<2)^z2)>>27;
        z2=((z2&4294967288U)<<2)^b;
        b=((z3<<13)^z3)>>21;
        z3=((z3&4294967280U)<<7)^b;
        b=((z4<<3)^z4)>>12;
        z4=((z4&4294967168U)<<13)^b;
        return (z1^z2^z3^z4);
    }
}

void srand(unsigned x)
{
    using namespace GenHelper;
    z1=x; 
    z2=(~x)^0x233333333U; 
    z3=x^0x1234598766U; 
    z4=(~x)+51;
}

int read()
{
    using namespace GenHelper;
    int a=rand_()&32767;
    int b=rand_()&32767;
    return a*32768+b;
}

const int maxn=2e7+1000;
const int INF=0x7fffffff;

int a[maxn],ls[maxn],rs[maxn],root;

inline void init(int n)
{
    static int stack[maxn];
    int top=0;
    a[0]=-INF;
    for (int i=1;i<=n;++i)
    {
        while (top && a[stack[top]]<=a[i]) ls[i]=stack[top--];
        rs[stack[top]]=i;
        stack[++top]=i;
    }
    root=stack[1];
}

inline uint64_t query(int l,int r)
{
    for (int x=root;;x=x<l?rs[x]:ls[x])
        if (x>=l && x<=r) return a[x];
}

int main()
{
    int n,m,s;
    scanf("%d%d%d",&n,&m,&s);
    srand(s);
    for (int i=1;i<=n;++i) a[i]=read();
    uint64_t ans=0;int l,r;
    init(n);
    for (int i=1;i<=m;++i)
    {
        l=read()%n+1;
        r=read()%n+1;
        ans+=l>r?query(r,l):query(l,r);
    }
    printf("%llu",ans);
}
```

## 树链剖分(LOJ模板题，带换根)

```cpp
#include <cstdio>
#include <cstring>
#include <algorithm>
#define ls l,m,o<<1
#define rs m+1,r,o<<1|1

typedef long long ll;

const int maxn=1e5+1000;

int val[maxn],top[maxn],fa[maxn],siz[maxn],son[maxn];
int f[maxn][30],depth[maxn],head[maxn],cnt,root,id[maxn],w[maxn];
ll sumv[maxn<<2],addv[maxn<<2],tot,n;

struct Edge
{
    int to,next;
}edge[maxn<<1];

inline void _add(int u,int v)
{
    edge[++cnt].next=head[u];
    edge[cnt].to=v;
    head[u]=cnt;
}

inline void add(int u,int v)
{
    _add(u,v);_add(v,u);
}

void dfs(int u,int fa,int dep)
{
    depth[u]=dep;
    ::fa[u]=fa;
    siz[u]=1;int maxs=-1;
    for (int i=head[u];i;i=edge[i].next)
    {
        int v=edge[i].to;
        if (v!=fa)
        {
            dfs(v,u,dep+1);
            siz[u]+=siz[v];
            if (siz[v]>maxs) son[u]=v,maxs=siz[v];
        }
    }
}

void dfs(int u,int topf)
{
    id[u]=++tot;
    w[tot]=val[u];
    top[u]=topf;
    if (!son[u]) return;
    dfs(son[u],topf);
    for (int i=head[u];i;i=edge[i].next)
    {
        int v=edge[i].to;
        if (v!=fa[u] && v!=son[u]) dfs(v,v);
    }
}

inline void pushup(int o){sumv[o]=sumv[o<<1]+sumv[o<<1|1];}

inline void pushdown(int o,int len)
{
    addv[o<<1]+=addv[o];
    addv[o<<1|1]+=addv[o];
    sumv[o<<1]+=addv[o]*(len-(len>>1));
    sumv[o<<1|1]+=addv[o]*(len>>1);
    addv[o]=0;
}

void build(int l,int r,int o)
{
    if (l==r) {sumv[o]=w[l];return;}
    int m=(l+r)>>1;
    build(ls);build(rs);
    pushup(o);
}

void update(int L,int R,int c,int l,int r,int o)
{
    if (L<=l && R>=r)
    {
        addv[o]+=c;
        sumv[o]+=c*(r-l+1);
        return;
    }
    pushdown(o,r-l+1);
    int m=(l+r)>>1;
    if (L<=m) update(L,R,c,ls);
    if (R> m) update(L,R,c,rs);
    pushup(o);
}

ll query(int L,int R,int l,int r,int o)
{
    if (L<=l && R>=r) return sumv[o];
    int m=(l+r)>>1;
    pushdown(o,r-l+1);
    ll tot=0;
    if (L<=m) tot+=query(L,R,ls);
    if (R> m) tot+=query(L,R,rs);
    return tot;
}

inline void AddRoute(int u,int v,int k)
{
    while (top[u]!=top[v])
    {
        if (depth[top[u]]<depth[top[v]]) u^=v^=u^=v;
        update(id[top[u]],id[u],k,1,n,1);
        u=fa[top[u]];
    }
    if (depth[u]>depth[v]) u^=v^=u^=v;
    update(id[u],id[v],k,1,n,1);
}

inline ll QueryRoute(int u,int v)
{
    ll ans=0;
    while (top[u]!=top[v])
    {
        if (depth[top[u]]<depth[top[v]]) u^=v^=u^=v;
        ans+=query(id[top[u]],id[u],1,n,1);
        u=fa[top[u]];
    }
    if (depth[u]>depth[v]) u^=v^=u^=v;
    ans+=query(id[u],id[v],1,n,1);
    return ans;
}

inline int lca(int u,int v)
{
    while (top[u]!=top[v])
    {
        if (depth[top[u]]<depth[top[v]]) u^=v^=u^=v;
        u=fa[top[u]];
    }
    return depth[u]>depth[v]?v:u;
}

inline int get_son(int u,int prec)
{
    for (int i=20;~i;--i)
        if (depth[f[u][i]]>depth[prec]) u=f[u][i];
    return u;
}

inline void AddSubTree(int u,int k)
{
    if (u==root) update(1,n,k,1,n,1);
    else if (lca(u,root)!=u)
        update(id[u],id[u]+siz[u]-1,k,1,n,1);
    else
    {
        update(1,n,k,1,n,1);
        int v=get_son(root,u);
        update(id[v],id[v]+siz[v]-1,-k,1,n,1);
    }
}

inline ll QuerySubTree(int u)
{
    if (u==root) return query(1,n,1,n,1);
    else if (lca(u,root)!=u)
        return query(id[u],id[u]+siz[u]-1,1,n,1);
    else
    {
        int v=get_son(root,u);
        return query(1,n,1,n,1)-query(id[v],id[v]+siz[v]-1,1,n,1);
    }
        // return query(1,n,1,n,1)-query(id[u]+1,id[u]+siz[u]-1,1,n,1);
        // // else return query(1,n,1,n,1);
}

inline void get_anc()
{
    for (int i=1;i<=n;++i) f[i][0]=fa[i];
    for (int i=1;i<=20;++i)
        for (int u=1;u<=n;++u)
            f[u][i]=f[f[u][i-1]][i-1];
}

int main()
{
    // freopen("tree2.in","r",stdin);
    // freopen("tree20.out","w",stdout);
    int m;
    scanf("%d",&n);root=1;
    for (int i=1;i<=n;++i) scanf("%d",val+i);
    for (int i=1,tmp;i<n;++i) scanf("%d",&tmp),add(tmp,i+1);
    dfs(1,0,1);dfs(1,1);build(1,n,1);get_anc();
    scanf("%d",&m);
    for (int i=1,opt,u,v,k;i<=m;++i)
    {
        scanf("%d",&opt);
        switch (opt)
        {
            case 1:scanf("%d",&root);break;
            case 2:scanf("%d%d%d",&u,&v,&k);AddRoute(u,v,k);break;
            case 3:scanf("%d%d",&u,&k);AddSubTree(u,k);break;
            case 4:scanf("%d%d",&u,&v);printf("%lld\n",QueryRoute(u,v));break;
            case 5:scanf("%d",&u);printf("%lld\n",QuerySubTree(u));break;
        }
    }
}
```



## LCT

```cpp
#include <cstdio>
#define ls(x) (ch[x][0])
#define rs(x) (ch[x][1])

const int maxn=3e5+1000;

int val[maxn],ch[maxn][2],fa[maxn],rev[maxn],s[maxn];

inline void swap(int &a,int &b){a^=b^=a^=b;}
inline bool nroot(int x){return ls(fa[x])==x || rs(fa[x])==x;}
inline void pushup(int x){s[x]=s[ls(x)]^s[rs(x)]^val[x];}

inline void pushr(int x)
{
    swap(ls(x),rs(x));rev[x]^=1;
}

inline void pushdown(int x)
{
    if (rev[x]) pushr(ls(x)),pushr(rs(x)),rev[x]=0;
}

inline void rotate(int x)
{
    int y=fa[x],z=fa[y],k=ch[y][1]==x;
    if (nroot(y)) ch[z][ch[z][1]==y]=x;
    fa[x]=z;
    ch[y][k]=ch[x][k^1];
    if (ch[x][k^1]) fa[ch[x][k^1]]=y;
    ch[x][k^1]=y;
    fa[y]=x;
    pushup(y);pushup(x);
}

inline void splay(int x)
{
    static int stack[maxn];
    int y=x,z=0;
    while (nroot(y)) stack[++z]=y,y=fa[y];
    stack[++z]=y;
    while (z) pushdown(stack[z--]);
    while (nroot(x))
    {
        y=fa[x],z=fa[y];
        if (nroot(y)) rotate(((ls(z)==y)^(ls(y)==x))?x:y);
        rotate(x);
    }
}

inline void access(int x)
{
    for (int y=0;x;y=x,x=fa[x])
        splay(x),ch[x][1]=y,pushup(x);
}

inline void makeroot(int x)
{
    access(x);splay(x);pushr(x);
}

inline int findroot(int x)
{
    access(x);
    splay(x);
    while (ls(x)) pushdown(x),x=ls(x);// 这里一定记得pushdown!
    splay(x);
    return x;
}

inline void link(int x,int y)
{
    makeroot(x);
    if (findroot(y)!=x) fa[x]=y;
}

inline void cut(int x,int y)
{
    makeroot(x);
    if (findroot(y)==x && fa[y]==x && !ch[y][0]) ch[x][1]=fa[y]=0,pushup(x);
}

inline void modify(int x,int y)
{
    splay(x);
    s[x]^=val[x];
    val[x]=y;
    s[x]^=val[x];
}

inline void split(int x,int y)
{
    makeroot(x);access(y);splay(y);
}

inline int query(int x,int y){split(x,y);return s[y];}

int main()
{
    int n,m;
    scanf("%d%d",&n,&m);
    for (int i=1;i<=n;++i)
        scanf("%d",val+i);
    for (int i=1,opt,x,y;i<=m;++i)
    {
        scanf("%d%d%d",&opt,&x,&y);
        switch (opt)
        {
            case 0:printf("%d\n",query(x,y));break;
            case 1:link(x,y);break;
            case 2:cut(x,y);break;
            case 3:modify(x,y);break;
        }
    }
}
```

## K-D Tree
```cpp
#include <cstdio>
#include <algorithm>

using std::nth_element;
using std::min;
using std::max;

const int maxn=1e6+100;
const int INF=0x3f3f3f3f;

struct Point
{
    int x,y;
    Point(){}
    Point(int x,int y):x(x),y(y){}
};

struct KDTree
{
    static constexpr double alpha=0.75;
    int root;
    struct Node
    {
        int val,siz,sum;
        int ch[2];
        Point mn,mx,now;
    }t[maxn];
    int top,has_rebuild,fa_rebuild,k_rebuild;
    KDTree(){t[0].mn=Point(INF,INF);t[0].mx=Point(-INF,-INF);}
    int trash[maxn],Trash;
    inline bool isbad(int o){return t[t[o].ch[0]].siz>t[o].siz*alpha || t[t[o].ch[1]].siz>t[o].siz*alpha;}
    inline void pushup(int o)
    {
        t[o].siz=t[t[o].ch[0]].siz+t[t[o].ch[1]].siz+1;
        t[o].sum=t[t[o].ch[0]].sum+t[t[o].ch[1]].sum+t[o].val;
        t[o].mn.x=min({t[o].now.x,t[t[o].ch[0]].mn.x,t[t[o].ch[1]].mn.x});
        t[o].mn.y=min({t[o].now.y,t[t[o].ch[0]].mn.y,t[t[o].ch[1]].mn.y});
        t[o].mx.y=max({t[o].now.y,t[t[o].ch[0]].mx.y,t[t[o].ch[1]].mx.y});
        t[o].mx.x=max({t[o].now.x,t[t[o].ch[0]].mx.x,t[t[o].ch[1]].mx.x});
    }
    inline void dfs(int o)
    {
        if (!o) return;
        if (t[o].ch[0]) dfs(t[o].ch[0]);
        trash[++Trash]=o;
        if (t[o].ch[1]) dfs(t[o].ch[1]);
    }
    inline int _rebuild(int l,int r,int k)
    {
        if (l>r) return 0;
        int mid=(l+r)>>1;
        int x=trash[mid];
        if (l==r) 
        {
            t[x].ch[0]=t[x].ch[1]=0;
            t[x].mn=t[x].mx=t[x].now;
            t[x].siz=1;
            t[x].sum=t[x].val;
            return x;
        }
        if (k==1) nth_element(trash+l,trash+mid+1,trash+r+1,[this](const int x,const int y) ->bool {return t[x].now.x<t[y].now.x;});
        else nth_element(trash+l,trash+mid+1,trash+r+1,[this](const int x,const int y) ->bool {return t[x].now.y<t[y].now.y;});
        x=trash[mid];
        t[x].ch[0]=_rebuild(l,mid-1,k^1);
        t[x].ch[1]=_rebuild(mid+1,r,k^1);
        pushup(x);
        return x;
    }
    inline int rebuild(int o,int k)
    {
        Trash=0;
        dfs(o);
        if (Trash) return _rebuild(1,Trash,k);
        return 0;
    }
    void _insert(const Point& p,const int x,int& o,int k)
    {
        if (!o)
        {
            o=++top;
            t[o].sum=t[o].val=x;
            t[o].mn=t[o].mx=t[o].now=p;
            t[o].siz=1;
            return;
        }
        if (k&1) _insert(p,x,t[o].ch[t[o].now.x<=p.x],0);
        else _insert(p,x,t[o].ch[t[o].now.y<=p.y],1);
        pushup(o);
        if (isbad(t[o].ch[0])) has_rebuild=t[o].ch[0],fa_rebuild=o,k_rebuild=k^1;
        else if (isbad(t[o].ch[1])) has_rebuild=t[o].ch[1],fa_rebuild=o,k_rebuild=k^1;
    }
    inline void insert(const Point& p,const int x)
    {
        _insert(p,x,root,0);
        if (isbad(root)) root=rebuild(root,0);
        else if (has_rebuild) t[fa_rebuild].ch[t[fa_rebuild].ch[1]==has_rebuild]=rebuild(has_rebuild,k_rebuild);
        has_rebuild=fa_rebuild=k_rebuild=0;
    }
    inline bool check_range(int o,const Point& l,const Point& r)
    {
        return t[o].mn.x>=l.x && t[o].mn.y>=l.y && t[o].mx.x<=r.x && t[o].mx.y<=r.y;
    }
    inline bool check_point(int o,const Point& l,const Point& r)
    {
        return (t[o].now.x>=l.x && t[o].now.x<=r.x) && (t[o].now.y>=l.y && t[o].now.y<=r.y);
    }
    inline bool check_have(int o,const Point& l,const Point& r)
    {
        return !((t[o].mx.x<l.x) || (t[o].mx.y<l.y) || (t[o].mn.x>r.x) || (t[o].mn.y>r.y));
    }
    inline int query(const Point& l,const Point& r,int o)
    {
        if (!o) return 0;
        if (check_range(o,l,r)) return t[o].sum;
        int ans=0;
        if (check_point(o,l,r)) ans+=t[o].val;
        if (check_have(t[o].ch[0],l,r)) ans+=query(l,r,t[o].ch[0]);
        if (check_have(t[o].ch[1],l,r)) ans+=query(l,r,t[o].ch[1]);
        return ans;
    }
}T;

int main()
{
    int n;
    scanf("%d",&n);
    int opt,x,lastans=0;
    Point a,b;
    while (scanf("%d",&opt) && opt!=3 && "STO LCH OTZ")
    {
        if (opt==1)
        {
            scanf("%d%d%d",&a.x,&a.y,&x);
            a.x^=lastans;a.y^=lastans;x^=lastans;
            T.insert(a,x);
        }
        else
        {
            scanf("%d%d%d%d",&a.x,&a.y,&b.x,&b.y);
            a.x^=lastans;a.y^=lastans;b.x^=lastans;b.y^=lastans;
            printf("%d\n",lastans=T.query(a,b,T.root));
        }
    }
}
```

## 平衡树

### Treap

```cpp
struct Treap
{
    struct Node
    {
    	int v,s,r,cnt;
    	Node* ch[2];
    	Node(int v,int s=1):s(s),cnt(s),r(rand()*rand()),v(v){ch[0]=ch[1]=0;}
    	inline void pushup()
    	{
    		s=cnt;
    		if (ch[0]) s+=ch[0]->s;
    		if (ch[1]) s+=ch[1]->s;
    	}
    	inline int cmp(int x)
    	{
    		return v==x?-1:v<x;
    	}
    };
	private:
		int __prec,__succ;
	public:
		Treap():__prec(0),__succ(0),root(0){}
		Node *root;
		inline void rotate(Node* &o,int d)
		{
			Node *k=o->ch[d^1];
			o->ch[d^1]=k->ch[d];
			k->ch[d]=o;
			o->pushup();
			k->pushup();
			o=k;
		}
	
		inline void _insert(Node* &o,int x,int t=1)
		{
			if (!o){o=new Node(x,t);return;}
			int d=o->cmp(x);
			if (d==-1) {o->cnt+=t;o->s+=t;return;}
			_insert(o->ch[d],x,t);
			if (o->ch[d]->r > o->r) rotate(o,d^1);
			o->pushup();
		}
		inline void insert(int x,int k=1){_insert(root,x,k);}
	
		inline void _remove(Node* &o,int x)
		{
			if (!o) return;
			int d=o->cmp(x);
			if (d==-1) 
			{
				if (o->cnt>1) {--o->cnt;--o->s;return;}
				if (!(o->ch[0])){Node* k=o;o=o->ch[1];delete k;return;}
				else if (!(o->ch[1])){Node* k=o;o=o->ch[0];delete k;return;}
				else
				{
					int d2=(o->ch[0]->r > o->ch[1]->r);
					rotate(o,d2);
					_remove(o->ch[d2],x);
				}
			}
			else _remove(o->ch[d],x);
			if (o) o->pushup();
		}
		inline void remove(int x){_remove(root,x);}
	
		inline int _kth(Node* o,int k)
		{
			if (!o || k<=0) return INT_MIN;
			int s=o->ch[0]?o->ch[0]->s:0;
			if (k>=s+1 && k<=s+o->cnt) return o->v;
			if (k<=s) return _kth(o->ch[0],k);
			return _kth(o->ch[1],k-s-o->cnt);
		}
		inline int kth(int k){return _kth(root,k);}
	
		inline int _rank(Node* o,int x)
		{
			if (!o) return 1;
			int s=o->ch[0]?o->ch[0]->s:0;
			if (o->v==x) return s+1;
			if (o->v<x) return s+o->cnt+_rank(o->ch[1],x);
			return _rank(o->ch[0],x);
		}
		inline int rank(int x){return _rank(root,x);}
		
		inline int _count(Node* o,int x)
		{
			if (!o) return 0;
			int d=o->cmp(x);
			if (~d) return _count(o->ch[d],x);
			return o->cnt;
		}
		inline int count(int x){return _count(root,x);}
	
		inline void _prec(Node* o,int x)
		{
			if (!o) return;
			if (o->v<x) __prec=max(__prec,o->v);
			if (o->v>=x) _prec(o->ch[0],x);
				else _prec(o->ch[1],x);
		}
		inline void _succ(Node* o,int x)
		{
			if (!o) return;
			if (o->v>x) __succ=min(__succ,o->v);
			if (o->v<=x) _succ(o->ch[1],x);
				else _succ(o->ch[0],x);
		}
		inline int prec(int x){__prec=INT_MIN+1;_prec(root,x);return __prec;}
		inline int succ(int x){__succ=INT_MAX  ;_succ(root,x);return __succ;}
}T;
```

### 替罪羊树
```cpp
#include <cstdio>
#include <vector>
#include <cassert>

using std::vector;

const int maxn=1e5+1000;

struct Scapegoat_Tree
{
    int root,has_rebuild,fa;
#if __cplusplus >= 201103L
    static constexpr double alpha=0.75;
#else 
    static const double alpha=0.75;
#endif
    struct Node
    {
        int siz,cnt,val;
        bool deleted;
        int ch[2];
    }t[maxn];
    inline bool isbad(int x){return t[t[x].ch[0]].cnt>t[x].cnt*alpha || t[t[x].ch[1]].cnt>t[x].cnt*alpha;}
    inline void maintain(int x){t[x].siz=t[t[x].ch[0]].siz+t[t[x].ch[1]].siz+!t[x].deleted;t[x].cnt=t[t[x].ch[0]].cnt+t[t[x].ch[1]].cnt+1;}
    int pool[maxn],top;
    vector<int> tmp;
    Scapegoat_Tree(){for (int i=1;i<maxn;++i) pool[++top]=i;}
    void dfs(int o)
    {
        if (t[o].ch[0]) dfs(t[o].ch[0]);
        if (!t[o].deleted) tmp.push_back(o);
        else pool[++top]=o;
        if (t[o].ch[1]) dfs(t[o].ch[1]);
    }
    int _rebuild(int l,int r)
    {
        if (l>r) return 0;
        int mid=(l+r)>>1,o=tmp[mid];
        if (l==r) t[o].siz=t[o].cnt=1;
        t[o].ch[0]=_rebuild(l,mid-1);
        t[o].ch[1]=_rebuild(mid+1,r);
        maintain(o);
        return o;
    }
    int rebuild(int o)
    {
        tmp.clear();
        dfs(o);
        return _rebuild(0,tmp.size()-1);
    }
    void _insert(int x,int &o)
    {
        if (!o)
        {
            o=pool[top--];
            t[o].siz=t[o].cnt=1;
            t[o].val=x;
            t[o].deleted=false;
            t[o].ch[0]=t[o].ch[1]=0;
            return;
        }
        ++t[o].siz;++t[o].cnt;
        int d=(t[o].val<=x);
        _insert(x,t[o].ch[d]);
        // if (isbad(o)) o=rebuild(o);//need repair?
        if (isbad(t[o].ch[0])) has_rebuild=t[o].ch[0],fa=o;
        else if (isbad(t[o].ch[1])) has_rebuild=t[o].ch[1],fa=o;
        // maintain(o);
    }
    void insert(int x) 
    {
        _insert(x,root);
        if (has_rebuild) t[fa].ch[t[fa].ch[1]==has_rebuild]=rebuild(has_rebuild);
        fa=has_rebuild=0;
    }
    void _remove_kth(int x,int o)
    {
        if (t[t[o].ch[0]].siz+1==x && !t[o].deleted)
        {
            t[o].deleted=true;
            t[o].siz--;
            return;
        }
        --t[o].siz;
        if (t[t[o].ch[0]].siz>=x) _remove_kth(x,t[o].ch[0]);
        else _remove_kth(x-t[t[o].ch[0]].siz-!t[o].deleted,t[o].ch[1]);
        // maintain(o);
    }
    void remove_kth(int x)
    {
        _remove_kth(x,root);
        if (t[root].cnt*alpha>t[root].siz) root=rebuild(root);
        fa=has_rebuild=0;
    }
    int rank(int x)
    {
        int ans=1;
        int o=root;
        while (o)
            if (t[o].val>=x) o=t[o].ch[0];
            else ans+=t[t[o].ch[0]].siz+!t[o].deleted,o=t[o].ch[1]; 
        return ans;
    }
    int kth(int k)
    {
        int o=root;
        while (1)
        {
            int s=t[t[o].ch[0]].siz;
            if (!t[o].deleted && s+1==k) return t[o].val;
            if (s>=k) o=t[o].ch[0];
            else k-=s+!t[o].deleted,o=t[o].ch[1];
        }
    }
    void print(int o)
    {
        if (!o) return;
        print(t[o].ch[0]);
        // printf("sizls=%d sizrs=%d ,siz=%d\n",t[t[o].ch[0]].siz,t[t[o].ch[1]].siz,t[o].siz);
        if (!t[o].deleted) printf("%d ",t[o].val);
        print(t[o].ch[1]);
    }
    /* data */
}T;

int main()
{
    int n;
    scanf("%d",&n);
    for (int i=1,opt,x;i<=n;++i)
    {
        scanf("%d%d",&opt,&x);
        switch (opt)
        {
            case 1:
                T.insert(x);
                break;
            case 2:
                T.remove_kth(T.rank(x));
                break;
            case 3:
                printf("%d\n",T.rank(x));
                break;
            case 4:
                printf("%d\n",T.kth(x));
                break;
            case 5:
                printf("%d\n",T.kth(T.rank(x)-1));
                break;
            case 6:
                printf("%d\n",T.kth(T.rank(x+1)));
                break;
        }
        // T.print(T.root);putchar('\n');
        // T.root=T.rebuild(T.root);
    }
}
```

### Splay

```cpp
struct Splay
{
    int root,cnt;
    struct Node
    {
        int val,size,ff,ch[2],cnt;
    }t[maxn];

    inline void pushup(int x)
    {
        t[x].size=t[x].cnt+t[t[x].ch[0]].size+t[t[x].ch[1]].size;
    }

    inline void rotate(int x)
    {
        int y=t[x].ff,z=t[y].ff;
        int k=(t[y].ch[1]==x);
        t[z].ch[t[z].ch[1]==y]=x;
        t[x].ff=z;
        t[y].ch[k]=t[x].ch[k^1];
        t[t[x].ch[k^1]].ff=y;
        t[y].ff=x;
        t[x].ch[k^1]=y;
        pushup(y);pushup(x);
    }

    inline void splay(int x,int goal)
    {
        while (t[x].ff!=goal)
        {
            int y=t[x].ff,z=t[y].ff;
            if (z!=goal) rotate((t[z].ch[1]==y)^(t[y].ch[1]==x)?x:y);
            rotate(x);
        }
        if (!goal) root=x;
    }

    inline void insert(int x)
    {
        int u=root,ff=0;
        while (t[u].val!=x && u)
        {
            ff=u;
            u=t[u].ch[t[u].val<x];
        }
        if (u) {++t[u].cnt;splay(u,0);return;}
        u=++cnt;
        if (ff) t[ff].ch[x>t[ff].val]=u;
        t[u].ff=ff;
        t[u].ch[0]=t[u].ch[1]=0;
        t[u].size=1;
        t[u].cnt=1;
        t[u].val=x;
        splay(u,0);
    }

    inline void find(int x)
    {
        int u=root;
        if (!u) return;
        while (t[u].ch[t[u].val<x] && t[u].val!=x) u=t[u].ch[t[u].val<x];
        splay(u,0);
    }

    inline int Next(int x,int type)
    {
        find(x);
        int u=root;
        if (t[u].val>x && type) return u;
        if (t[u].val<x && !type) return u;
        u=t[u].ch[type];
        while (t[u].ch[type^1])u=t[u].ch[type^1];
        splay(u,0);
        return u;
    }

    inline void Delete(int x)
    {
        int prev=Next(x,0),succ=Next(x,1);
        splay(prev,0);splay(succ,prev);
        int u=t[succ].ch[0];
        if (t[u].cnt>1){--t[u].cnt;splay(u,0);return;}
            else t[succ].ch[0]=0;
    }

    inline int kth(int k)
    {
        int u=root;
        if (t[u].size<k) return inf;
        while (19260817)
        {
            int s=t[t[u].ch[0]].size;
            if (s+t[u].cnt<k)
                k-=s+t[u].cnt,u=t[u].ch[1];
            else if (k<=s)
                u=t[u].ch[0];
            else {splay(u,0);return t[u].val;}
        } 
    }
}T;
```

### 权值线段树

```cpp
void insert(int x,int l=1,int r=n,int o=1)
{
    ++sumv[o];
    if (l==r) return;
    int m=(l+r)>>1;
    if (x<=m) insert(x,ls);
        else insert(x,rs);
}

inline void remove(int x,int l=1,int r=n,int o=1)
{
    --sumv[o];
    if (l==r) return;
    int m=(l+r)>>1;
    if (x<=m) remove(x,ls);
        else remove(x,rs);
} 

inline int count(int x,int l=1,int r=n,int o=1)
{
    if (l==r) return sumv[o];
    int m=(l+r)>>1;
    if (x<=m) return count(x,ls);
        else return count(x,rs);
}

inline int rank(int x,int l=1,int r=n,int o=1)
{
    if (l==r) return 1;
    int m=(l+r)>>1;
    if (x<=m) return rank(x,ls);
        else return rank(x,rs)+sumv[o<<1];
}

inline int kth(int k,int l=1,int r=n,int o=1)
{
    if (l==r) return l;
    int m=(l+r)>>1;
    if (k<=sumv[o<<1]) return kth(k,ls);
        return kth(k-sumv[o<<1],rs);
}

inline int prec(int x)
{
    return kth(rank(x)-1);
}

inline int succ(int x)
{
    return kth(rank(x)+count(x));
}

inline int get_rnk(int i)
{
    return lower_bound(b+1,b+n+1,i)-b;
}

```



## 树套树

```cpp
struct Node
{
	int v,s,r,cnt;
	Node* ch[2];
	Node(int v,int s=1):s(s),cnt(s),r(rand()*rand()),v(v){ch[0]=ch[1]=0;}
	inline void pushup()
	{
		s=cnt;
		if (ch[0]) s+=ch[0]->s;
		if (ch[1]) s+=ch[1]->s;
	}
	inline int cmp(int x)
	{
		return v==x?-1:v<x;
	}
};

struct Treap
{
	private:
		int __prec,__succ;
	public:
		Treap():__prec(0),__succ(0),root(0){}
		Node *root;
		inline void rotate(Node* &o,int d)
		{
			Node *k=o->ch[d^1];
			o->ch[d^1]=k->ch[d];
			k->ch[d]=o;
			o->pushup();
			k->pushup();
			o=k;
		}
	
		inline void _insert(Node* &o,int x,int t=1)
		{
			if (!o){o=new Node(x,t);return;}
			int d=o->cmp(x);
			if (d==-1) {o->cnt+=t;o->s+=t;return;}
			_insert(o->ch[d],x,t);
			if (o->ch[d]->r > o->r) rotate(o,d^1);
			o->pushup();
		}
		inline void insert(int x,int k=1){_insert(root,x,k);}
	
		inline void _remove(Node* &o,int x)
		{
			if (!o) return;
			int d=o->cmp(x);
			if (d==-1) 
			{
				if (o->cnt>1) {--o->cnt;--o->s;return;}
				if (!(o->ch[0])){Node* k=o;o=o->ch[1];delete k;return;}
				else if (!(o->ch[1])){Node* k=o;o=o->ch[0];delete k;return;}
				else
				{
					int d2=(o->ch[0]->r > o->ch[1]->r);
					rotate(o,d2);
					_remove(o->ch[d2],x);
				}
			}
			else _remove(o->ch[d],x);
			if (o) o->pushup();
		}
		inline void remove(int x){_remove(root,x);}
	
		inline int _kth(Node* o,int k)
		{
			if (!o || k<=0) return INT_MIN;
			int s=o->ch[0]?o->ch[0]->s:0;
			if (k>=s+1 && k<=s+o->cnt) return o->v;
			if (k<=s) return _kth(o->ch[0],k);
			return _kth(o->ch[1],k-s-o->cnt);
		}
		inline int kth(int k){return _kth(root,k);}
	
		inline int _rank(Node* o,int x)
		{
			if (!o) return 1;
			int s=o->ch[0]?o->ch[0]->s:0;
			if (o->v==x) return s+1;
			if (o->v<x) return s+o->cnt+_rank(o->ch[1],x);
			return _rank(o->ch[0],x);
		}
		inline int rank(int x){return _rank(root,x);}
		
		inline int _count(Node* o,int x)
		{
			if (!o) return 0;
			int d=o->cmp(x);
			if (~d) return _count(o->ch[d],x);
			return o->cnt;
		}
		inline int count(int x){return _count(root,x);}
	
		inline void _prec(Node* o,int x)
		{
			if (!o) return;
			if (o->v<x) __prec=max(__prec,o->v);
			if (o->v>=x) _prec(o->ch[0],x);
				else _prec(o->ch[1],x);
		}
		inline void _succ(Node* o,int x)
		{
			if (!o) return;
			if (o->v>x) __succ=min(__succ,o->v);
			if (o->v<=x) _succ(o->ch[1],x);
				else _succ(o->ch[0],x);
		}
		inline int prec(int x){__prec=INT_MIN+1;_prec(root,x);return __prec;}
		inline int succ(int x){__succ=INT_MAX  ;_succ(root,x);return __succ;}
}tree[maxn<<2];

inline void insert_tree(Treap& t,Node* rt)
{
	if (!rt) return;
	t.insert(rt->v,rt->cnt);
	insert_tree(t,rt->ch[0]);
	insert_tree(t,rt->ch[1]);
}

inline void pushup(int o)
{
	insert_tree(tree[o],tree[o<<1].root);
	insert_tree(tree[o],tree[o<<1|1].root);
}

inline void build(int l,int r,int o)
{
	if (l==r) {tree[o].insert(a[l]);return;}
	int m=(l+r)>>1;
	build(ls);
	build(rs);
	pushup(o);
}

inline int rank(int L,int R,int k,int l,int r,int o)
{
	if (L<=l && R>=r) return tree[o].rank(k);
	const int m=(l+r)>>1;
	int tot=0;
	if (L<=m) tot+=rank(L,R,k,ls);
	if (R> m) tot+=rank(L,R,k,rs);
	if (L<=m && R>m) --tot;
	return tot;
}

inline int count(int L,int R,int x,int l,int r,int o)
{
	if (L<=l && R>=r) return tree[o].count(x);
	const int m=(l+r)>>1;
	int tot=0;
	if (L<=m) tot+=count(L,R,x,ls);
	if (R> m) tot+=count(L,R,x,rs);
	return tot;
}

inline void update(int p,int x,int l,int r,int o)
{
	tree[o].remove(a[p]);
	tree[o].insert(x);
	if (l==r) return;
	int m=(l+r)>>1;
	if (p<=m) update(p,x,ls);
		else update(p,x,rs);
}

inline int prec(int L,int R,int x,int l,int r,int o)
{
	if (L<=l && R>=r) return tree[o].prec(x);
	int m=(l+r)>>1,ans=INT_MIN+1;
	if (L<=m) ans=prec(L,R,x,ls);
	if (R>m)  ans=max(ans,prec(L,R,x,rs));
	return ans;
}

inline int succ(int L,int R,int x,int l,int r,int o)
{
	if (L<=l && R>=r) return tree[o].succ(x);
	int m=(l+r)>>1,ans=INT_MAX;
	if (L<=m) ans=succ(L,R,x,ls);
	if (R> m) ans=min(ans,succ(L,R,x,rs));
	return ans;
}

inline int kth(int L,int R,int k,int n)
{
	int l=1,r=tree[1].root->s;
	while (l<=r)
	{
		int m=(l+r)>>1;
		int t=tree[1].kth(m);
		int K=rank(L,R,t,1,n,1);
		int cnt=count(L,R,t,1,n,1);
		if (K==k && !cnt) return succ(L,R,t,1,n,1);
		if (K+cnt-1>=k && K<=k && cnt) return t;
		if (cnt)
			if (K+cnt-1<k) l=m+1;
				else r=m-1;
		else
			if (K<k) l=m+1;
				else r=m-1;
	}
	return INT_MAX;
}
```

## DSU ON TREE（CF600E）

```cpp
#include <cstdio>
#include <cstring>
#include <algorithm>

const int maxn=1e5+1000;

typedef long long ll;

ll ans[maxn],tot;
int head[maxn],_cnt,n,mx;
int fa[maxn],siz[maxn],son[maxn],cnt[maxn],col[maxn],isson[maxn];

struct Edge
{
    int to,next;
}edge[maxn<<1];

inline int max(int a,int b){return a<b?b:a;}

inline void add(int u,int v)
{
    edge[++_cnt].next=head[u];
    edge[_cnt].to=v;
    head[u]=_cnt;
}

void dfs(int u,int f)
{
    fa[u]=f;
    siz[u]=1;
    int maxs=-1;
    for (int i=head[u];i;i=edge[i].next)
    {
        int v=edge[i].to;
        if (v!=f)
        {
            dfs(v,u);siz[u]+=siz[v];
            if (siz[v]>maxs) maxs=siz[v],son[u]=v;
        } 
    }
    if (son[u]) isson[son[u]]=true;
}

void AddAns(int u,int x,int Son)
{
    if((cnt[col[u]]+=x)>mx) mx=cnt[col[u]],tot=col[u];
        else if(cnt[col[u]]==mx) tot+=(ll)col[u];
    for (int i=head[u];i;i=edge[i].next)
    {
        int v=edge[i].to;
        if (v!=fa[u] && v!=Son) AddAns(v,x,Son);
    }
}

void solve(int u,int f,int keep)
{
    for (int i=head[u];i;i=edge[i].next)
    {
        int v=edge[i].to;
        if (v!=f && v!=son[u]) solve(v,u,0);
    }
    if (son[u]) solve(son[u],u,1);
    AddAns(u,1,son[u]);
    ans[u]=tot;
    // int t=0;
    // for (int i=1;i<=n;++i) t=max(t,cnt[i]);
    // for (int i=1;i<=n;++i) if (cnt[i]==t) ans[u]=ans[u]+(ll)i;
    if (!keep) AddAns(u,-1,0),mx=tot=0;
}

int main()
{
    scanf("%d",&n);
    for (int i=1;i<=n;++i) scanf("%d",col+i);
    for (int i=1,u,v;i<n;++i)
    {
        scanf("%d%d",&u,&v);
        add(u,v);add(v,u);
    }
    dfs(1,1);
    solve(1,0,0);
    for (int i=1;i<=n;++i) printf("%lld ",ans[i]);
}
```



## 珂朵莉树

```cpp
struct Node
{
    int l,r;
    mutable int x;
    Node(int l,int r=-1,int x=0):l(l),r(r),x(x){}
    bool operator< (const Node &nd) const {return l<nd.l;}
};

typedef long long ll;
typedef set<Node>::iterator It;

struct Cmp
{
    bool operator() (const It& a,const It& b) const {return a->x<b->x;}
};

int a[maxn];
int seed,n,m,vmax;
set<Node> s;

inline void swap(int& a,int& b){a^=b^=a^=b;}

inline int Rand()
{
    int ret=seed;
    seed=((ll)seed*7+13)%1000000007;
    return ret;
}

inline void init()
{
    for (int i=1;i<=n;++i)
        s.insert(Node(i,i,a[i]));
    s.insert(Node(n+1,n+1,0));
}

inline ll pow_mod(ll a,ll b,ll p)
{
    ll ans=1%p;a%=p;
    for (;b;b>>=1)
    {
        if (b&1) ans=ans*a%p;
        a=a*a%p;
    }
    return ans;
}

inline It split(int pos)
{
    It it=s.lower_bound(Node(pos)); 
    if (it!=s.end() && it->l==pos) return it;
    --it;
    int l=it->l,r=it->r,v=it->x;
    s.erase(it);
    s.insert(Node(l,pos-1,v));
    return s.insert(Node(pos,r,v)).first;
}

inline void assign(int l,int r,int x)
{
    It it2=split(r+1),it1=split(l);
    s.erase(it1,it2);
    s.insert(Node(l,r,x));
}

inline void add(int l,int r,int x)
{
    It it2=split(r+1),it1=split(l);
    for (It i=it1;i!=it2;++i) i->x+=x;
}

inline int kth(int l,int r,int x)
{
    vector<It> v;
    It it2=split(r+1),it1=split(l);
    for (It i=it1;i!=it2;++i) v.push_back(i);
    sort(v.begin(),v.end(),Cmp());
    for (int i=0;i<v.size();++i)
    {
        // assert(v[i]->r>=v[i]->l);//Assertion Failed...Fixed....
        x-=v[i]->r-v[i]->l+1;
        if (x<=0) return v[i]->x;
    }
    return -1;
}

inline int power(int l,int r,int x,int p)
{
    ll ans=0;
    It it2=split(r+1),it1=split(l);
    for (It it=it1;it!=it2;++it) ans=((ll)ans+(ll)pow_mod(it->x,x,p)*(it->r-it->l+1)%p)%p;
    return ans;
}
```

