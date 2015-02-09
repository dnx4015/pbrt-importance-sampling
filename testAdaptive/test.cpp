/*
ID: diana.n1
PROG: 
LANG: C++
*/

#include <cstdio>
#include <cstring>
#include <cmath>
#include <string>
#include <algorithm>
#include <sstream>
#include <cstdlib>
using namespace std;

#define IOR(x) freopen(x,"r",stdin);
#define IOW(x) freopen(x,"w",stdout);
#define DEBUG if(0)

#define i64 long long
#define u64 unsigned long long
#define eps 1e-10

#define REP(i,n) for(int i=0;i<n;i++)
#define FOR(i,a,b) for(int i=a;i<=b;i++)
#define FORIT(it,p) for(__typeof(p.end()) it=p.begin();it!=p.end();it++)
#define INC(i,a,b, inc) for(int i=a;i<=b;i+=inc)

#define mset(p,v) memset(p,v,sizeof(p))
#define all(x) x.begin(), x.end()
#define mp make_pair
#define fst first
#define snd second
#define pb push_back

#define SSTR( x ) dynamic_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()

using namespace std;

inline float calcDistance(std::pair<int, float>p1, 
						  std::pair<int, float>p2, 
				          std::pair<int, float>p){
	float dist = fabs((p2.fst - p1.fst)*(p1.snd - p.snd)
					- (p1.fst - p.fst)*(p2.snd - p1.snd));
	return dist / sqrt((p2.fst - p1.fst)*(p2.fst - p1.fst)
			+ (p2.snd - p1.snd)*(p2.snd - p1.snd));
}

int main(){
    pair<int, float>p1 = mp(0,0);
    pair<int, float>p2 = mp(32,1000);
    pair<int, float>p = mp(1,31.249);
	float dist = calcDistance(p1, p2, p);
    printf("%f\n",dist);
    return 0;
}

