/* This software is copyright (C) 2009-2012 Malcolm Sharpe, Levent Tuncel, and
 * Tor Myklebust.  All rights reserved.
 *
 * You may modify and use this software for non-commercial purposes.  You may
 * redistribute the source code of modified versions of this software provided
 * that:
 *  (a) Your redistributions of the source code contains the above copyright
 *      notice, this list of conditions, and the following disclaimer, and
 *  (b) The modified software is prominently marked as modified.
 *
 * THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESSED OR IMPLIED WARRANTIES,
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
 * FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL MALCOLM
 * SHARPE, LEVENT TUNCEL, OR TOR MYKLEBUST BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
 * OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/* This software has some modifications made by Fernandes. Ferreira, Francoa and Schouery.
 * It was only modified in terms of input/output. 
 * None of the original heuristics implemented by Sharp, Tuncel and Myklebust where modified.
 */

#include <algorithm>
#include <assert.h>
#include <cctype>
#include <execinfo.h>
#include <cfloat>
#include <bitset>
#include <sys/types.h>
#include <sys/resource.h>
#include <sys/signal.h>
#include <unistd.h>
#include <signal.h>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <queue>
#include <set>
#include <utility>
#include <vector>
#include <string.h>
using namespace std;

#define FR(i,a,b) for(int i=(a);i<(b);++i)
#define FOR(i,n) FR(i,0,n)
#define RF(i,a,b) for(int i=(a)-1;i>=(b);--i)
#define ROF(i,n) RF(i,n,0)
#define ROFI(i,v) ROF(i,(int)(v).size())
#define BEND(v) (v).begin(),(v).end()
#define MP make_pair
#define PB push_back
#define CLR(x,a) memset(x,a,sizeof(x));
#define FORALL(i,v) for(typeof((v).end())i=(v).begin();i!=(v).end();++i)
#define dump(x) cerr << #x << " = " << (x) << endl
#ifdef assert
#undef assert
#endif
#define assert(x) \
    do { \
      bool foo = x; \
      if (!foo) { \
        printf("Assertion failed (%s:%i): %s\n", __FILE__, __LINE__, #x); \
        dump_backtrace(); \
        exit(-1); \
      } \
    } while (0)
          

long long get_cpu_usecs() {
  rusage r;
  getrusage(RUSAGE_SELF, &r);
  return            r.ru_utime.tv_usec
      + 1000000LL * r.ru_utime.tv_sec;
}

void dump_backtrace() {
  static void *p[1024];
  static char buf[1024];
  int k = backtrace(p,1000);
  sprintf(buf, "/usr/bin/addr2line -f -C -e /proc/%i/exe", getpid());
  FILE *f = popen(buf, "w");
  FOR(i,k) fprintf(f, "%p\n", p[i]);
  fflush(f);
  fclose(f);
}

const double inf = 1e20;

// make this fairly big, since the input never uses very many decimal places
const double eps = 1e-3;

#define MAXN 44444
#define MAXM 44444+1

template <typename T>
struct matrix {
  int sx, sy;
  T *data;

  matrix() : sx(0), sy(0), data(0) {}

  matrix(int sx, int sy) : sx(sx), sy(sy) { data = new T[sx*sy]; }
  matrix(const matrix &m) : sx(m.sx), sy(m.sy) {
    data = new T[sx*sy];
    FOR(i, sx*sy) data[i] = m.data[i];
  }
  matrix(matrix &&m) : sx(m.sx), sy(m.sy), data(m.data) {
    m.data = 0;
    m.sx = m.sy = 0;
  }
  ~matrix() { delete[] data; }
  matrix &operator=(const matrix &m) {
    if (this != &m) {
      delete[] data;
      sx = m.sx; sy = m.sy;
      data = new T[sx*sy];
      FOR(i, sx*sy) data[i] = m.data[i];
    }
    return *this;
  }

  T *operator[](size_t k) { return data + k * sy; }
  size_t size() { return sx; }
};

int n,m;
matrix<double> R;
double N[MAXN];
double delta;
vector<double> pi;
vector<int> asst;

// options
bool do_good = 0; // repeatedly do pc() || ust() || crazy() then C and Pi.
bool do_maxr = 0; // do maxr; otherwise assign everyone to null product
bool do_guru_et_al = 0; // do my distortion of the Guruswami et al. log-approx algo
bool do_gurutun = 0; // do Levent's enhancement of the Guruswami et al. log-approx algo
bool do_dk = 0; // do dobson-kalish
bool do_torsdk = 0; // do Tor's implementation of Dobson-Kalish
bool do_pidk = 0; // do price-space dobson-kalish
bool do_pd = 0; // do price-down heuristic
bool do_pu = 0; // do price-up heuristic
bool do_ufc = 0; // do update-from-choices heuristic
bool do_ust = 0; // do update subtrees heuristic
bool do_nldk = 0; // do non-local dobson-kalish
bool do_crazy = 0; // do crazy random search
bool do_only_small_steps = 0; // reduce pi by minimum to change asst
bool do_prefer_small_steps = 0; // prefer small steps in pricechange
bool do_recompute_asst = 0;
bool do_recompute_pi = 0;
bool do_iterate_operator = 0; // repeat recomputation while it helps
bool do_ub2 = 0; // use upper bound #2 in the searches
bool do_expensive_first = 0; // in BnB, first try assigning to high-priced items
bool do_cheap_first = 0; // in BnB, first try assigning to low-priced items
bool do_poor_first = 0; // in BnB, first try assigning consumers with low maxR
bool do_rich_first = 0; // in BnB, first try assigning consumers with high maxR
bool do_random_order = 0; // in BnB, try options in random order
bool do_dump_soln = 0; // output the optimal assignments
bool do_dump_dir = 0; // output the results of searching along various directions
bool do_dump_pi = 0; // output the optimal prices
bool do_trace = 0; // trace the execution of the searches
bool paranoid = 0; // verify frequently (slow for DK)
bool verbose = 0; // print a lot more
double bnb_tol = 0.05; // percentage within optimal for branch-and-bound

bool torsdk();

void parse_args(int argc, char *argv[]) {
  FR(i,1,argc) {
    char *arg = argv[i];

    if (0==strcmp(arg,"--good")) do_good = 1;
    if (0==strcmp(arg,"--maxr")) do_maxr = 1;
    if (0==strcmp(arg,"--guru-et-al")) do_guru_et_al = 1;
    if (0==strcmp(arg,"--gurutun")) do_gurutun = 1;
    if (0==strcmp(arg,"--dk")) do_dk = 1;
    if (0==strcmp(arg,"--torsdk")) do_torsdk = 1;
    if (0==strcmp(arg,"--pidk")) do_pidk = 1;
    if (0==strcmp(arg,"--pd")) do_pd = 1;
    if (0==strcmp(arg,"--pu")) do_pu = 1;
    if (0==strcmp(arg,"--pc")) do_pd = do_pu = 1;
    if (0==strcmp(arg,"--ufc")) do_ufc = 1;
    if (0==strcmp(arg,"--ust")) do_ust = 1;
    if (0==strcmp(arg,"--nldk")) do_nldk = 1;
    if (0==strcmp(arg,"--crazy")) do_crazy = 1;
    if (0==strcmp(arg,"--only-small-steps")) do_only_small_steps = 1;
    if (0==strcmp(arg,"--prefer-small-steps")) do_prefer_small_steps = 1;
    if (0==strcmp(arg,"--recompute-asst")) do_recompute_asst = 1;
    if (0==strcmp(arg,"--recompute-pi")) do_recompute_pi = 1;
    if (0==strcmp(arg,"--iterate-operator")) do_iterate_operator = 1;
    if (0==strcmp(arg,"--ub2")) do_ub2 = 1;
    if (0==strcmp(arg,"--expensive-first")) do_expensive_first = 1;
    if (0==strcmp(arg,"--cheap-first")) do_cheap_first = 1;
    if (0==strcmp(arg,"--poor-first")) do_poor_first = 1;
    if (0==strcmp(arg,"--rich-first")) do_rich_first = 1;
    if (0==strcmp(arg,"--random-order")) do_random_order = 1;
    if (0==strcmp(arg,"--dump-soln")) do_dump_soln = 1;
    if (0==strcmp(arg,"--dump-dir")) do_dump_dir = 1;
    if (0==strcmp(arg,"--dump-pi")) do_dump_pi = 1;
    if (0==strcmp(arg,"--trace")) do_trace = 1;
    if (0==strcmp(arg,"--paranoid")) paranoid = 1;
    if (0==strcmp(arg,"--verbose")) verbose = 1;

    if (0==strcmp(arg,"--help")) {
      printf("Usage: %s [options]\n"
	     "\n"
	     "Input is taken from stdin and written to stdout.\n"
	     "\n"
	     "Valid options:\n"
	     "  --maxr: initialize with MAXR heuristic\n"
	     "  --guru-et-al: do my distortion of the Guruswami et al. log-approx algo\n"
	     "  --gurutun: do Levent's enhancement of the Guruswami et al. log-approx algo\n"
	     "  --dk: use Dobson-Kalish heuristic\n"
	     "  --torsdk: use alternate implementation of Dobson-Kalish heuristic\n"
	     "  --pidk: use price-space Dobson-Kalish heuristic\n"
	     "  --pd: use pricedown heuristic\n"
	     "  --pu: use priceup heuristic\n"
	     "  --pc: use pricechange heuristic (same as --pd and --pu)\n"
	     "  --ufc: use update-from-choices heuristic\n"
	     "  --ust: use update subtrees heuristic\n"
	     "  --nldk: use non-local Dobson-Kalish heuristic\n"
	     "  --crazy: do crazy random search\n"
	     "  --only-small-steps: use local price reductions in pricedown\n"
	     "  --prefer-small-steps: prefer small steps in pricechange\n"
	     "  --recompute-asst: recompute optimal assignment\n"
	     "  --recompute-pi: recompute optimal prices\n"
	     "  --iterate-operator: repeatedly perform recomputation\n"
	     "  --dump-soln: output the optimal assignment\n"
	     "  --dump-dir: output the second-best result for each coordinate search direction\n"
	     "  --dump-pi: output the optimal prices\n"
	     "  --trace: trace the execution of the searches\n"
	     "  --paranoid: verify frequently (slower)\n"
	     "  --verbose: print a lot more\n"
	     ,argv[0]);
      exit(0);
    }
  }
}

// true if heuristics should be completely silent
//   (used for branch-and-bound)
bool silent=0;

// get consumer i's surplus if assigned to product j
inline double getsurp(int i, int j) {
  return R[i][j] - pi[j];
}

// total consumer surplus
double gettotsurp() {
  double ret = 0;
  FOR(i,n) ret += N[i] * getsurp(i, asst[i]);
  return ret;
}

// total welfare
double getwelfare() {
  double ret = 0;
  FOR(i,n) ret += N[i] * R[i][asst[i]];
  return ret;
}

// objective value
double getobj() {
  double obj = 0;
  FOR(i,n) {
    //printf("asst[%d] = %d (price %d)\n",i,asst[i],pi[asst[i]]);
    obj += N[i] * pi[asst[i]];
  }
  return obj;
}

// get the number of products that are sold to someone
bool issold[MAXM];
int getnsold() {
  CLR(issold,0);
  issold[0] = 1;
  FOR(i,n) issold[asst[i]] = 1;
  
  int ret = 0;
  FOR(j,m) ret += issold[j];
  return ret;
}

// check that the current prices and assignment are feasible
void verify() {
  assert(pi[0] == 0);

  FOR(i,n) {
    FOR(j,m) if (asst[i] != j) {
      if (getsurp(i,asst[i])+eps < getsurp(i,j)) {
	printf("  consumer %d wants to switch from %d to %d\n",i,asst[i],j);
      }
      assert(getsurp(i, asst[i])+eps >= getsurp(i, j));
    }
  }

  FOR(j,m) {
    if (pi[j] > inf) printf("%f\n", pi[j]);
    assert(pi[j] <= inf+eps);
  }
}

// compute an optimal assignment given the current prices
void recompute_asst() {
  FOR(i,n) {
    int bestj = 0;
    double bestsurp = 0;

    FOR(j,m) {
      if (getsurp(i,j) > bestsurp+eps
	|| (fabs(getsurp(i,j)-bestsurp)<eps && pi[j] > pi[bestj]+eps)) {
	bestj = j;
	bestsurp = getsurp(i,j);
      }
    }

    asst[i] = bestj;
  }
}

// compute arc costs for the price computation
// if a node j0 is specified, compute only the arc costs for the in-arcs of j0
bitset<MAXM> tedg[MAXM];
matrix<double> adj;

void build_graph(int j0=-1) {
  FOR(j,m) if (j0==-1 || j==j0) FOR(k,m) adj[k][j] = inf;

  FOR(i,n) {
    int j = asst[i];

    // this consumer is unassigned (different from being assigned to null)
    if (j==-1) continue;

    if (j0==-1 || j==j0) FOR(k,m) {
      // reduce cost by pi so that all arc costs are nonnegative
      double rc = (R[i][j] - R[i][k]) - (pi[j] - pi[k]);
      assert(rc >= -eps);
      adj[k][j] = min(adj[k][j], rc);
    }
  }
}

// compute shortest path distances on unmarked nodes,
// assuming marked nodes have distance zero
bool mark[MAXM];
double y[MAXM];
void dijkstra() {
  FOR(j,m) y[j] = mark[j] ? 0 : inf;

  // start with path distances coming from marked nodes
  FOR(k,m) if (!mark[k]) {
    FOR(j,m) if (mark[j]) {
      y[k] = min(y[k], y[j] + adj[j][k]);
    }
  }

  while (1) {
    int j = -1;
    FOR(k,m) if (!mark[k]) if (j==-1 || y[k] < y[j]) j = k;
    if (j==-1) break;
    mark[j] = 1;

    if (y[j]+eps >= inf) break; // all remaining nodes have infinite distance

    FOR(k,m) {
      y[k] = min(y[k], y[j] + adj[j][k]);
    }
  }
}

// compute optimal prices given the current assignment
void recompute_pi() {
  // compute shortest path distances reduced by
  // the previous prices
  build_graph();
  //verify();

  CLR(mark,0);
  mark[0] = 1;
  dijkstra();

  // update pi from the reduced shortest path distances
  FOR(j,m) {
    assert(y[j] >= -eps);
    if (y[j]+eps >= inf) pi[j] = inf;
    else pi[j] += y[j];
  }

  verify();
}

// the event of a customer switching to a product
struct Evt {
  int i; // the customer
  int j; // the product the customer switches to
  double pj; // the price at which the consumer switches into the product
  double palt; // the price of the consumer's best alternative to j

  Evt(int i, int j, double pj, double palt) : i(i), j(j), pj(pj), palt(palt) {}
};

// events are ordered as they occur as a price is raised
bool operator<(const Evt &a, const Evt &b) {
  // if (fabs((a.pj-pi[a.j]) - (b.pj-pi[b.j])) > eps) return a.pj+pi[b.j] < b.pj+pi[a.j];

  // a.j and b.j must agree on whether they initially are infinite,
  // or else the search can't be done properly
  bool ajinf = pi[a.j]+eps >= inf, bjinf = pi[b.j]+eps >= inf;

  if (ajinf) {
    // in the infinite case, the prices are considered to have "started at
    // the same infinite point"
    assert(bjinf);

    if (fabs(a.pj - b.pj) > eps) return a.pj < b.pj;
  } else {
    // in the non-infinite case, the prices are moving at the same rate
    assert(!bjinf);

    if (fabs((a.pj-pi[a.j]) - (b.pj-pi[b.j])) > eps) return a.pj+pi[b.j] < b.pj+pi[a.j];
  }

  return a.pj+b.palt < b.pj+a.palt;
}

// the choice of consumer i buying product j, in order of preference
// given current prices
struct Choice {
  int i;
  int j;
  Choice(int i, int j) : i(i), j(j) {}
};

// can compare choices of a single consumer
int cmpchoice(int i, int j, int k) {
  if (fabs(getsurp(i,j)-getsurp(i,k)) > eps) {
    if (getsurp(i,j) < getsurp(i,k)) return -1;
    return 1;
  }
  if (fabs(pi[j] - pi[k]) > eps) {
    if (pi[j] < pi[k]) return -1;
    return 1;
  }
  return 0;
}
bool operator<(const Choice &a, const Choice &b) {
  assert(a.i == b.i);
  return cmpchoice(a.i, a.j, b.j) < 0;
}

// the preferences of each consumer at current prices
vector<Choice> choices[MAXN];
int getalt(int i, int j) {
  if (j != choices[i][m-1].j) return choices[i][m-1].j;
  assert(m >= 2);
  return choices[i][m-2].j;
}

// compute, for each consumer, its preference order at current prices
void compute_choices() {
  FOR(i,n) choices[i].clear();
  FOR(i,n) {
    FOR(j,m) {
      choices[i].PB(Choice(i,j));
    }
    sort(BEND(choices[i]));
  }
}

// try to change a single price if it helps and return true
// otherwise return false
double pcdir[MAXM];
bool update() {
  int bestj = -1;
  double bestpj = 0;
  double oldobj = getobj();
  double curobj = getobj(), bestobj = curobj;
  double bestratio = 0;

  compute_choices();

  // try changing prices for all non-null products
  FR(j,1,m) {
    pcdir[j] = 0;

    // indifference events for all customers, between product j and the best alternative
    vector<Evt> indiff;

    FOR(i,n) {
      int k = getalt(i,j);
      double price = R[i][j] - (R[i][k] - pi[k]);

      indiff.PB(Evt(i,j,price,pi[k]));
    }

    sort(BEND(indiff));

    // gradually lower the price of j from infinity,
    // computing the resulting objective value
    double baseobj = 0;

    FOR(i,n) {
      baseobj += N[i] * pi[getalt(i,j)];
    }

    int nevt = indiff.size();
    double Mj = 0;
    ROF(k,nevt+1) {
      double pj = inf;

      if (k < nevt) {
	int i = indiff[k].i;
	pj = indiff[k].pj;
	double palt = indiff[k].palt;

	assert(fabs(palt-pi[getalt(i,j)]) < eps);

	Mj += N[i];
	baseobj -= palt * N[i];
      }

      if (!do_pu && pj > pi[j]+eps) continue;
      if (!do_pd && pj < pi[j]-eps) continue;

      double newobj = baseobj + Mj * pj;

      if (fabs(newobj-oldobj) > eps) pcdir[j] = max(pcdir[j], newobj);

      bool better;
      if (do_prefer_small_steps) {
	better = newobj > oldobj+eps
	    && (newobj-oldobj)/(1+fabs(pj-pi[j])) > bestratio;
      } else {
	better = newobj > bestobj+eps;
      }

      if (better) {
	bestobj = newobj;
	bestj = j;
	bestpj = pj;
	bestratio = (newobj-oldobj)/(1+fabs(pj-pi[j]));
      }

      if (pj+eps < pi[j] && do_only_small_steps) break;
    }
  }

  if (bestj == -1) return 0;

  int j = bestj;
  if (!silent) printf("  changed price of %d from %lf to %lf\n", j, pi[j], bestpj);
  pi[j] = bestpj;

  // fix assignment
  FOR(i,n) {
    double altsurp = getsurp(i, getalt(i,j)),
	surpj = getsurp(i, j);
    int newasst;
    if (altsurp+eps < surpj || (fabs(altsurp-surpj)<eps && pi[j] > pi[getalt(i,j)]+eps)) {
      newasst = j;
    } else {
      newasst = getalt(i,j);
    }
    if (asst[i] != newasst && !silent && verbose) printf("  reassigned consumer %d\n", i);
    asst[i] = newasst;
  }

  if (!silent) printf("    (giving obj = %lf; expected = %lf)\n", getobj(), bestobj);
  verify();
  if (!do_prefer_small_steps) assert(fabs(bestobj-getobj()) < eps);

  return 1;
}

// try to change several prices taken from the top of
// a consumer's preference list
int invchoice[MAXM];
//int wantpfx[MAXN][MAXM+1], wantsfx[MAXN][MAXM+1];
vector<vector<int> > wantpfx, wantsfx;

bool update_from_choices() {
  int besti = -1;
  int bestc = -1;
  double bestdp = 0;
  double curobj = getobj(), bestobj = curobj;

  compute_choices();

  // try changing prices for the most-liked products
  // of those consumers not assigned to null
  FOR(i,n) if (asst[i]) {
    // determine the ranks of the products
    FOR(c,m) {
      invchoice[choices[i][c].j] = c;
    }

    // precompute the most preferred products for each consumer
    // from a prefix and suffix of i's choice list
    //CLR(wantpfx,-1);
    //CLR(wantsfx,-1);
    wantpfx = vector<vector<int> >(n, vector<int>(m+1, -1));
    wantsfx = vector<vector<int> >(n, vector<int>(m+1, -1));

    FOR(ii,n) {
      // first compute the highest rank in the choice list that
      // can be achieved
      FOR(d,m) {
	int k = choices[ii][d].j;
	wantsfx[ii][invchoice[k]] = d;
	wantpfx[ii][invchoice[k]+1] = d;
      }

      FOR(d,m) {
	wantpfx[ii][d+1] = max(wantpfx[ii][d+1], wantpfx[ii][d]);
	wantsfx[ii][m-d-1] = max(wantsfx[ii][m-d-1], wantsfx[ii][m-d]);
      }

      // now replace the ranks with products
      FOR(d,m+1) {
	if (wantpfx[ii][d] != -1) wantpfx[ii][d] = choices[ii][wantpfx[ii][d]].j;
	if (wantsfx[ii][d] != -1) wantsfx[ii][d] = choices[ii][wantsfx[ii][d]].j;
      }
    }

    FOR(c,m) {
      // we can't change the price of the null product, and we
      // can't handle changing the prices of unbought products
      // (note that unbought products will be at the bottom of every
      // consumer's preferences)
      if (choices[i][m-c-1].j == 0 || pi[choices[i][m-c-1].j]+eps >= inf) break;

      // indifference events for all consumers, between the top c+1 products
      // and the best alternative
      vector<Evt> indiff;

      double baseobj = 0;
      FOR(ii,n) {
	int besttop = wantsfx[ii][m-c-1];
	int bestalt = wantpfx[ii][m-c-1];

	assert(besttop != -1);
	assert(bestalt != -1);

	int j = besttop, k = bestalt;

	double price = R[ii][j] - (R[ii][k] - pi[k]);

	indiff.PB(Evt(ii,j,price,pi[k]));

	baseobj += N[ii] * pi[k];
      }

      sort(BEND(indiff));

      // gradually lower the prices from infinity,
      // computing the resulting objective value
      int nevt = indiff.size();
      double Mtop = 0;
      ROF(k,nevt+1) {
	double dp = inf;

	if (k < nevt) {
	  int ii = indiff[k].i;
	  dp = indiff[k].pj - pi[indiff[k].j];
	  double palt = indiff[k].palt;

	  Mtop += N[ii];
	  baseobj -= palt * N[ii];
	  baseobj += pi[indiff[k].j] * N[ii];
	}

	double newobj = baseobj + Mtop * dp;

	if (newobj > bestobj+eps) {
	  besti = i;
	  bestobj = newobj;
	  bestc = c;
	  bestdp = dp;
	}
      }
    }
  }

  if (besti == -1) return 0;

  int i = besti;
  int c = bestc;
  printf("  changed price of consumer %d's top %d choices by %lf\n",
    i, c+1, bestdp);

  FOR(d,c+1) {
    pi[choices[i][m-d-1].j] += bestdp;
  }

  // fix assignment
  recompute_asst();

  printf("    (giving obj = %lf; expected = %lf)\n", getobj(), bestobj);
  verify();
  assert(fabs(bestobj-getobj()) < eps);

  return 1;
}

//// SUBTREE SEARCH DIRECTIONS
// try changing the prices of entire subtrees of the shortest path tree
void dfs(int);
bool line_search(double *);

bool update_subtrees() {
  recompute_pi();
  build_graph();

  vector<int> old_asst = asst;
  vector<double> old_pi = pi;
  double old_obj = getobj();
  vector<int> best_asst = asst;
  vector<double> best_pi = pi;
  double best_obj = getobj();
  int best_j = -1;

  for (int j = 1; j < m; j++) if (pi[j] < inf) {
    CLR(mark,0);
    dfs(j);
    vector<double> dir(m);
    FOR(k,m) dir[k] = mark[k];
    asst = old_asst; pi = old_pi;
    if (line_search(&dir[0])) {
      if (getobj() > best_obj) {
        best_asst = asst;
        best_pi = pi;
        best_obj = getobj();
        best_j = j;
      }
    }
  }
  if (best_obj > old_obj) {
    pi = best_pi;
    asst = best_asst;
    printf("  changed price of subtree rooted at %d\n", best_j);
    printf("    (giving obj = %lf)\n", getobj());
    return 1;
  } else {
    printf("ust failed\n");
    return 0;
  }
}

// mark as many non-descendants of a node in the shortest path tree as can occur
// for any such tree
void dfs(int j) {
  if (mark[j]) return;
  mark[j] = 1;
  FOR(k,m) if (fabs(adj[j][k]) < eps) dfs(k);
}
void marknondesc(int j) {
  CLR(mark,0);
  mark[j] = 1;
  dfs(0);
  FOR(k,m) if (pi[k]-eps >= inf) mark[k] = 1;
  mark[j] = 0;
}

//// NON-LOCAL VARIANT OF DOBSON-KALISH

// which products can have their prices vary
bool varying[MAXM];
// which products are allowed to lose customers
bool highend[MAXM];

vector<double> halt_time;
vector<int> nldk_asst;
vector<double> nldk_M;
double nldk_pi(int j) {
  assert(!varying[j]);

  return pi[j] + halt_time[j];
}

// reassignment events
struct NLDKEvt {
  int i; // the consumer whose choice is causing pi[j] to be blocked
  int j; // the product being blocked
  int k; // the product that i is indifferent to compared to product j
  double t; // the time (in unit price) when the event occurs

  // the price of j when the event occurs
  double pj() const { return pi[j] + t; }

  // the price of k when the event occurs
  double pk() const { return nldk_pi(k); }
};
bool operator<(const NLDKEvt &a, const NLDKEvt &b) {
  if (fabs(a.t-b.t) > eps) return a.t < b.t;

  // prefer to process reassignment events that help first
  return a.pj() - a.pk() < b.pj() - b.pk();
}
bool operator>(const NLDKEvt &a, const NLDKEvt &b) {
  return b < a;
}

// use a min-heap
typedef priority_queue<NLDKEvt, vector<NLDKEvt>, greater<NLDKEvt> > nldk_evt_queue;

nldk_evt_queue nldk_evts;
void nldk_scan(int k) {
  assert(!varying[k]);
  if (nldk_pi(k)+eps >= inf) return;
  FOR(i,n) {
    int j = nldk_asst[i];
    double t = (R[i][j] - pi[j]) - (R[i][k] - nldk_pi(k));
    if (highend[j]) {
      NLDKEvt evt;
      evt.i = i;
      evt.j = j;
      evt.k = k;
      evt.t = t;

      nldk_evts.push(evt);
    } else {
      halt_time[j] = min(halt_time[j], t);
    }
  }
}

// do non-local dobson-kalish (up)
double nldk_curobj, nldk_bestobj;
int nldk_bestj;
double nldk_bestpj;
double nldk_nvary, nldk_newobj;
double nldk_t;
void nldk_init(int j) {
  if (verbose) fprintf(stderr, "  RUNNING NLDK FOR PRODUCT %d\n", j);

  nldk_evts = nldk_evt_queue();
  halt_time = vector<double>(m,inf);
  nldk_asst = asst;

  nldk_M = vector<double>(m,0);
  FOR(i,n) nldk_M[asst[i]] += N[i];

  // allow all prices for non-null products
  // vary initially if they can,
  // and select a single product to possibly lose customers.
  // this generalizes 1988 Dobson-Kalish.
  CLR(highend,0);
  highend[j] = 1;

  marknondesc(j);

  FOR(k,m) varying[k] = !mark[k];

  nldk_nvary = 0;
  FOR(k,m) if (varying[k]) nldk_nvary += nldk_M[k];

  FOR(k,m) if (!varying[k]) {
    halt_time[k] = 0;
    nldk_scan(k);
  }

  FOR(k,m) assert(pi[k]+eps < inf || !varying[k]);

  nldk_newobj = nldk_curobj;
}
// if target is set, terminate when that objective is achieved
bool nldk_exec_event(int j, double target=-1) {
  while (1) {
    int k = -1;
    FOR(l,m) if (varying[l] && !highend[l]) {
      if (k==-1 || halt_time[l] < halt_time[k]) {
	k = l;
      }
    }

    // discard invalid events
    while (nldk_evts.size()) {
      NLDKEvt evt = nldk_evts.top();

      assert(varying[evt.j]);
      if (nldk_asst[evt.i] == evt.j) {
	break;
      }

      nldk_evts.pop();
    }

    if (k != -1 && (nldk_evts.size()==0 || halt_time[k]-eps <= nldk_evts.top().t)) {
      if (verbose) fprintf(stderr, "    HALTING %d\n", k);

      nldk_t = halt_time[k];

      nldk_newobj += nldk_t * nldk_M[k];
      varying[k] = 0;
      nldk_nvary -= nldk_M[k];

      nldk_scan(k);
    } else if (nldk_evts.size()) {
      NLDKEvt evt = nldk_evts.top(); nldk_evts.pop();

      nldk_t = evt.t;

      double theobj = nldk_newobj + nldk_t*nldk_nvary;
      if (theobj > nldk_bestobj+eps) {
	nldk_bestobj = theobj;
	nldk_bestj = j;
	nldk_bestpj = evt.pj();
      }

      // abort if we've reached the best objective value we found before
      if (fabs(theobj-target) < eps) {
	return 0;
      }

      if (verbose) fprintf(stderr, "    REASSIGNING %d from %d to %d\n", evt.i,evt.j,evt.k);

      assert(highend[evt.j]);

      nldk_asst[evt.i] = evt.k;

      nldk_M[evt.j] -= N[evt.i];
      nldk_nvary -= N[evt.i];

      nldk_M[evt.k] += N[evt.i];

      nldk_newobj -= N[evt.i] * pi[evt.j];
      nldk_newobj += N[evt.i] * evt.pk();
    } else {
      if (verbose) fprintf(stderr, "    DONE RAISING PRICE OF %d\n", j);

      nldk_t = inf;
      assert(fabs(nldk_M[j]) < eps);
      assert(fabs(nldk_nvary) < eps);

      if (nldk_newobj > nldk_bestobj+eps) {
	nldk_bestobj = nldk_newobj;
	nldk_bestj = j;
	nldk_bestpj = inf;
      }
      return 0;
    }

    assert(nldk_nvary >= -eps);

    return 1;
  }
}
bool nldk() {
  nldk_curobj = nldk_bestobj = getobj();
  nldk_bestj=-1;
  nldk_bestpj=0;

  FR(j,1,m) if (pi[j]+eps < inf) {
    nldk_init(j);
    while (nldk_exec_event(j));
  }

  if (nldk_bestj==-1) return 0;
  assert(nldk_bestobj > getobj());

  nldk_init(nldk_bestj);
  while (nldk_exec_event(nldk_bestj, nldk_bestobj));

  FOR(k,m) {
    if (varying[k]) {
      if (nldk_t+eps < inf) {
	pi[k] += nldk_t;
      } else {
	pi[k] = nldk_t;
      }
    } else {
      if (pi[k]+eps >= inf) assert(halt_time[k] == 0);
      pi[k] += halt_time[k];
    }
  }
  asst = nldk_asst;
  if (paranoid) verify();

  if (verbose) fprintf(stderr, "  EXPECT OBJECTIVE %lf\n", nldk_bestobj);
  if (verbose) fprintf(stderr, "  RECEIVE OBJECTIVE %lf\n", getobj());

  assert(fabs(nldk_bestobj-getobj()) < eps);
  return 1;
}

//// DOBSON-KALISH (1988)
// do dobson-kalish
bool dk() {
  recompute_pi();
  build_graph();

  double curobj = getobj(), bestobj = curobj;
  int besti=-1, bestk=-1;

  FOR(i,n) if (asst[i]) {
    int j = asst[i];
    FOR(k,m) if (k != j) {
      if (fabs(adj[k][j]) < eps && fabs(getsurp(i,j)-getsurp(i,k)) < eps) {
	// the arc kj lies on a shortest path from 0 to j
	// and this is the consumer that makes it happen,
	// so try reassigning consumer i to product k
	if (verbose) printf("  trying to reassign consumer %d from %d to %d: ",i,j,k);

	marknondesc(j);
	
	asst[i] = k;
	if (paranoid) verify();

	// ensure only j and its descendants are unmarked,
	// and compute their shortest path distances
	build_graph(j);
	dijkstra();
	assert(fabs(y[k]) < eps);

	if (verbose) printf("(kj: %lf) ",adj[k][j]);

	// compute the new objective value
	double newobj = 0;
	FOR(ii,n) {
	  int jj = asst[ii];
	  newobj += N[ii] * (pi[jj]+y[jj]);
	}
	
	if (verbose) printf("(%lf -> %lf) ",pi[j],pi[j]+y[j]);
	if (verbose) printf("%lf\n",newobj);

	if (newobj > bestobj) {
	  bestobj = newobj;
	  besti = i;
	  bestk = k;
	}

	if (paranoid) {
	  // check that we computed the new prices and objective correctly
	  vector<double> oldy;
	  FOR(jj,m) oldy.PB(y[jj]);
	  vector<double> oldpi = pi;
	  recompute_pi();
	  
	  FOR(jj,m) {
	    if (verbose) printf("  EXPECT %lf+%lf=%lf RECEIVED %lf\n",
	      oldpi[jj],y[jj],oldpi[jj]+y[jj],pi[jj]);
	    assert((y[jj]+eps >= inf && pi[jj]+eps >= inf) ||
	      fabs(oldpi[jj]+y[jj]-pi[jj]) < eps);
	  }
	  if (verbose) printf("%lf = %lf\n",newobj,getobj());
	  assert(fabs(newobj-getobj()) < eps);
	  pi = oldpi;
	}

	// clean up
	asst[i] = j;
	build_graph(j);
	if (paranoid) build_graph(); // clean up after the paranoia above
	if (paranoid) verify();
      }
    }
  }

  if (besti == -1) return 0;

  if (!silent) printf("  reassigned consumer %d from %d to %d\n", besti, asst[besti], bestk);
  asst[besti] = bestk;
  recompute_pi();

  if (!silent) printf("    (giving obj = %lf)\n", getobj());
  assert(fabs(bestobj-getobj()) < eps);

  return 1;
}

// do price-space Dobson-Kalish
//   i.e. try raising each price by epsilon
const double pidk_eps = 0.1;
bool pidk() {
  double bestobj = getobj();
  int bestj=-1;

  FR(j,1,m) {
    vector<int> oldasst = asst;
    vector<double> oldpi = pi;

    pi[j] += pidk_eps;

    recompute_asst();
    recompute_pi();

    double theobj = getobj();
    if (theobj > bestobj+eps) {
      bestobj = theobj;
      bestj = j;
    }

    asst = oldasst;
    pi = oldpi;
  }

  if (bestj==-1) return 0;

  pi[bestj] += pidk_eps;

  recompute_asst();
  recompute_pi();

  assert(fabs(bestobj-getobj()) < eps);

  return 1;
}

//// RANDOM SEARCH
// just try random stuff
bool oldcrazy() {
  int ncrazy = 100*m;

  double curobj = getobj();

  double maxp = 0;

  FOR(i,n) FOR(j,m) maxp = max(maxp, R[i][j]);

  FOR(a,ncrazy) {
    vector<int> oldasst = asst;
    vector<double> oldpi = pi;

    FOR(j,m) {
      if (pi[j] >= inf-eps) {
	if (rand() <= RAND_MAX/m) pi[j] = rand() * maxp / RAND_MAX;
      } else {
	pi[j] *= rand() * 2.0 / RAND_MAX;
      }
    }

    recompute_asst();
    recompute_pi();

    if (getobj() > curobj+eps) {
      return 1;
    }

    asst = oldasst;
    pi = oldpi;
  }

  return 0;
}

//#define VERBOSE_LINE_SEARCH
bool line_search(double *dir) {
    #ifdef VERBOSE_LINE_SEARCH
    printf("beginning line_search with dir:\n");
    FOR(j, m) printf("%6.2f\t", dir[j]);
    printf("\n");
    printf("reservation prices:\n");
    FOR(i, n) {
      printf("<%4i> ", asst[i]);
      FOR(j,m) printf("%6.2f\t", R[i][j]);
      printf("\n");
    }
    printf("current prices:\n");
    FOR(j, m) printf("%6.2f\t", pi[j]);
    printf("\n");
    printf("\n");
    #endif

    double curobj = getobj();

    double mint = -1e99, maxt = 1e99;
    FOR(j,m) if (dir[j] < 0) maxt = min(maxt, -pi[j] / dir[j]);
    FOR(j,m) if (dir[j] > 0) mint = max(mint, -pi[j] / dir[j]);

    vector<pair<double, int> > prods;
    FOR(j, m) prods.push_back(make_pair(-dir[j], j));
    sort(prods.begin(), prods.end());

    struct break_event {
      double alpha;
      double objchange;
      double deltachange;
      int cust, prod, fromprod;

      break_event(double alpha, int cust, int from, int to, double *dir)
          : alpha(alpha) {
        double f = getsurp(cust, from) - alpha*dir[from];
        double t = getsurp(cust, to) - alpha*dir[to];
        objchange = t-f;
        deltachange = dir[to] - dir[from];
        this->cust = cust;
        this->prod = to;
        this->fromprod = from;
      }

      bool operator<(const break_event &b) const {
        if (alpha - b.alpha) return alpha < b.alpha;
        if (objchange - b.objchange) return objchange > b.objchange;
        if (deltachange - b.deltachange) return deltachange < b.deltachange;
        if (cust - b.cust) return cust < b.cust;
        return prod < b.prod;
      }
    };

    vector<break_event> events;

    vector<int> hull[n];
    FOR(i, n) {
      // Carefully do Graham's scan to find all of the products bought by i on
      // the line described by dir.
      int best = 0;
      int jj = 1;
      for (jj = 1; prods[jj].first == prods[0].first; jj++)
        if (getsurp(i,prods[jj].second) > getsurp(i, prods[best].second)
         || (getsurp(i,prods[jj].second) == getsurp(i,prods[best].second) &&
             pi[prods[jj].second] > pi[prods[best].second]))
          best = jj;
      hull[i].push_back(prods[best].second);
      while (jj < (int)prods.size()) {
        int jjj = jj+1;
        int best = jj;
        while (jjj < (int)prods.size() && prods[jjj].first == prods[jj].first) {
          if (getsurp(i, prods[jjj].second) > getsurp(i, prods[best].second)
           || (getsurp(i,prods[jjj].second) == getsurp(i,prods[best].second) &&
               pi[prods[jjj].second] > pi[prods[best].second]))
            best = jjj;
          jjj++;
        }
        jj = best;
        while (hull[i].size() >= 2u) {
          int j1 = hull[i][hull[i].size()-2], j2 = hull[i][hull[i].size()-1],
              j3 = prods[jj].second;
          double x1 = dir[j1], x2 = dir[j2], x3 = dir[j3];
          double y1 = getsurp(i,j1), y2 = getsurp(i,j2), y3 = getsurp(i,j3);
          double cross = (y3-y2)*(x2-x1) - (y2-y1)*(x3-x2);
          if (cross > 0) break;
          hull[i].pop_back();
        }
        hull[i].push_back(prods[jj].second);
        jj = jjj;
      }
      FOR(jj, hull[i].size()-1) {
        double m1 = dir[hull[i][jj]], m2 = dir[hull[i][jj+1]];
        double b1 = getsurp(i,hull[i][jj]), b2 = getsurp(i,hull[i][jj+1]);
        double alpha = (b2-b1)/(m2-m1);
        events.push_back(break_event(alpha, i, hull[i][jj], hull[i][jj+1], dir));
#ifdef VERBOSE_LINE_SEARCH
        {
          int hijj = hull[i][jj];
          int hijjj = hull[i][jj+1];
          printf("cust %i breaks from %i(%f+%ft) to %i(%f+%ft) at time %f\n",
  	    i,
  	    hijj, R[i][hijj]-pi[hijj], -dir[hijj],
  	    hijjj, R[i][hijjj]-pi[hijjj], -dir[hijjj],
  	    alpha);
        }
#endif
      }
    }
#ifdef VERBOSE_LINE_SEARCH
    FOR(i,n) {
      printf("hull %i(%f+%ft): ", i, pi[i],dir[i]);
      FOR(j,hull[i].size()) printf(" %i", hull[i][j]);
      printf("\n");
    }
#endif
    
    sort(events.begin(), events.end());

#ifdef VERBOSE_LINE_SEARCH
    FOR(i, events.size()) {
      printf("event: %i from %i to %i alpha %f bc %f mc %f\n",
        events[i].cust, events[i].fromprod, events[i].prod,
        events[i].alpha, events[i].objchange, events[i].deltachange);
    }
#endif

    double coef[n], b[n];
#ifdef VERBOSE_LINE_SEARCH
    int prodno[n];
    FOR(i,n) prodno[i] = hull[i][0];
#endif
    double slope = 0, yint = 0;
    FOR(i,n) coef[i] = dir[hull[i][0]], b[i] = pi[hull[i][0]];
    FOR(i,n) slope += N[i] * coef[i], yint += N[i] * b[i];
    double best = -1.0/0.0, bestt = 0;

    int hitmint = 0;
    FOR(z, events.size()) {
      break_event &e = events[z];
      double alpha = e.alpha;
#ifdef VERBOSE_LINE_SEARCH
      if (alpha >= mint && alpha <= maxt) {
        printf("considering point alpha = %f:\n", alpha);
        FOR(j,m) printf("pi[%i] = %f\n", j, pi[j] + alpha*dir[j]);
        printf("welfares:\n");
        FOR(i,n) {
          double welf[m];
          FOR(j,m) {
            double fff = R[i][j] - pi[j] - alpha*dir[j];
            welf[j] = fff;
            if (fff < -999) fff = -999;
            if (fff > 9999) fff = 9999;
            if (j == prodno[i]) printf("<%4.0f>", fff);
            else printf(" %4.0f ", fff);
          }
          printf("\n");
          FOR(j,m) if (welf[j] > welf[prodno[i]] + 1e-2) {
            printf("BUGGY %i %i\n", i, j);
          }
        }
        double profit = 0, welfare = 0;
        FOR(i,n) profit += N[i] * (pi[prodno[i]] + alpha*dir[prodno[i]]);
        FOR(i,n) welfare += N[i] * (R[i][prodno[i]] - pi[prodno[i]] - alpha*dir[prodno[i]]);
        printf("postulated profit is %f\n", slope * alpha + yint);
        printf("profit is %f; total welfare is %f\n", profit, welfare);
      }
#endif

#define TRY(x) \
        if (fabs(x) > 1e-4 && slope * (x) + yint > best) \
          best = slope * (x) + yint, bestt = (x);
      if (alpha > maxt) {
        TRY(maxt);
        break;
      }
      if (alpha >= mint) {
        if (!hitmint) TRY(mint);
        hitmint = true;
        TRY(alpha);
      }
      slope -= N[e.cust] * coef[e.cust];
      yint -= N[e.cust] * b[e.cust];
      slope += N[e.cust] * (coef[e.cust] = dir[e.prod]);
      yint += N[e.cust] * (b[e.cust] = pi[e.prod]);
#ifdef VERBOSE_LINE_SEARCH
      prodno[e.cust] = e.prod;
#endif
      if (alpha >= mint) TRY(alpha);
#undef TRY
    }

    if (best > curobj + 1e-4) {
      FOR(j,m) pi[j] += bestt * dir[j];
      recompute_asst();
      return true;
    }
    return false;
}

bool crazy() {
  long long before = get_cpu_usecs();

  vector<vector<bool> > adj(m, vector<bool>(m, false));
  FOR(i,n) {
    int j = asst[i];
    if (j) FOR(k, m) if (k && j != k && R[i][j]-pi[j]==R[i][k]-pi[k]) {
      adj[j][k] = 1;
    }
  }
  FOR(j,m) FOR(k,m) if (adj[j][k]) {
    double dir[m];
    FOR(a, m) dir[a] = (j==a) - (k==a);
    if (line_search(dir)) return true;
    break;
  }

  long long after = get_cpu_usecs();
  printf("spent %lli usecs in crazy()\n", after-before);

  return false;
}

// collection of heuristics
void do_heuristic_init() {
  // initialize all non-null product prices to infinity
  pi = vector<double>(m, inf);
  pi[0] = 0; // null product always has price zero

  // initially assign all customers to the null product
  asst = vector<int>(n,0);

  if (do_maxr) {
    // if using maxr heuristic, initially assign each
    // consumer to the product which he/she likes the most,
    // then compute prices
    printf("initializing with MAXR heuristic\n");

    pi = vector<double>(m, 0);
    FOR(i,n) {
      FOR(j,m) {
	if (R[i][j] > R[i][asst[i]]) {
	  asst[i] = j;
	}
      }
    }
    recompute_pi();
  }

  if (do_guru_et_al) {
    // this modified version of Guruswami et al. log-approx algo
    // sets all prices to the same value, which is given by each
    // consumer's maximum reservation price
    printf("initializing with Guruswami et al.\n");

    double bestp = 0, bestobj = 0;

    FOR(i,n) {
      double maxr = 0;
      FOR(j,m) maxr = max(maxr, R[i][j]);

      pi = vector<double>(m, maxr);
      pi[0] = 0;
      recompute_asst();
      if (do_recompute_pi) recompute_pi();
      
      if (getobj() > bestobj+eps) {
	bestobj = getobj();
	bestp = maxr;
      }
    }

    pi = vector<double>(m, bestp);
    pi[0] = 0;
    recompute_asst();
    if (do_recompute_pi) recompute_pi();
  } else if (do_gurutun) {
    // do Guruswami et al.'s log-approx algo with Levent's enhancement

    // initially, we did not worry about lower bounds on prices,
    // assuming that any violated lower bound would be reflected in decreased objective

    // while this is true, this may cause us to miss good moves
    // so now we check the lower bounds
    printf("initializing with Levent's enhanced Guruswami et al.\n");

    vector<int> prod;
    FR(j,1,m) prod.PB(j);

    pi = vector<double>(m, inf);
    pi[0] = 0;
    recompute_asst();

    while (1) {
      double bestp = inf, bestobj = 0;
      double prevobj = getobj();

      vector<double> pilb(m,0);

      FOR(i,n) if (asst[i]) {
	FOR(j,m) {
	  pilb[j] = max(pilb[j], R[i][j] - getsurp(i,asst[i]));
	}
      }

      FOR(i,n) if (asst[i] == 0) {
	double maxr = 0;
	FORALL(j,prod) maxr = max(maxr, R[i][*j]);

	vector<int> oldasst = asst;
	vector<double> oldpi = pi;
	FR(j,1,m) if (maxr >= pilb[j]-eps) pi[j] = maxr;
	recompute_asst();

	if (getobj() > bestobj+eps) {
	  bestobj = getobj();
	  bestp = maxr;
	}
	pi = oldpi;
	asst = oldasst;
      }

      if (bestobj <= prevobj+eps) break;

      FR(j,1,m) if (bestp >= pilb[j]-eps) pi[j] = bestp;
      recompute_asst();

      prod.clear();

      set<int> bought;
      FOR(i,n) bought.insert(asst[i]);
      FR(j,1,m) if (!bought.count(j)) prod.PB(j);

      printf("  MOVE TO OBJECTIVE %lf (nprod = %d)\n", getobj(), (int)prod.size());
    }
  }

  if (paranoid) verify();
}


void do_heuristics(int *count=0) {
  if (paranoid) verify();

  while (1) {
    if (!silent && count) printf("count = %d, obj = %lf, welfare = %lf\n", *count, getobj(), getwelfare());

    // try heuristics and continue if one works; otherwise stop
    bool change = 0;
    if (do_pd || do_pu) if (update()) change = 1;
    if (do_ufc) if (update_from_choices()) change = 1;
    if (do_ust) if (update_subtrees()) change = 1;
    if (do_dk) if (dk()) change = 1;
    if (do_torsdk) if (torsdk()) change = 1;
    if (do_pidk) if (pidk()) change = 1;
    if (do_nldk) if (nldk()) change = 1;
    if (do_crazy) if (crazy()) change = 1;
    if (!change) break;

    double oldobj;
    do {
      if (!silent) printf("  recomputing... ");
      oldobj = getobj();

      // recompute optimal assignment
      if (do_recompute_asst) recompute_asst();

      // recompute optimal prices
      if (do_recompute_pi) recompute_pi();

      if (!silent) printf("to obj = %lf\n", getobj());
    } while (do_iterate_operator && oldobj+eps < getobj());

    if (count) ++*count;
  }
}

void do_good_heuristics(int *count=0) {
  while (update_subtrees());
  do_pd = do_pu = true;
  while (1) {
    if (!silent && count)
      printf("count = %d, obj = %lf, welfare = %lf\n",
             *count, getobj(), getwelfare());
    bool change = false;
    while (1) {
      bool ch = torsdk();
      if (!ch) break;
      change |= ch;
      if (count) ++*count;
    }
    if (change) recompute_pi();
    else change = update() || update_subtrees() || crazy();
    if (!change) break;
    recompute_asst();
    recompute_pi();
    if (count) ++*count;
  }
}

void dumpsoln() {
  FOR(i,n) {
    int j = asst[i];
    printf("  Consumer %d buys %d for %lf (R = %lf)\n", i, j, pi[j], R[i][j]);
  }
}

void dumppi() {
  FOR(j,m) {
    printf("  Product %d costs %lf\n", j, pi[j]);
  }
}

void dumpdir() {
  FR(j,1,m) {
    printf("  Direction e_%d gives objective %lf\n", j, pcdir[j]);
  }
}

long long begin_time;

void last_words() {
  long long fin_time = get_cpu_usecs();
  printf("\n");
  printf("Products sold: %d / %d\n", getnsold()-1, m-1);
  printf("Total welfare: %lf\n", getwelfare());
  printf("Total consumer surplus: %lf\n", gettotsurp());
  assert(fabs(gettotsurp() + getobj() - getwelfare()) < eps);
  printf("Final objective value: %lf\n", getobj());
  printf("cpu usecs: %lli\n", fin_time - begin_time);
}

void sigdie(int sig) {
  printf("Dying due to signal %i\n", sig);
  last_words();
  exit(0);
}

void setup_signal_handlers() {
  signal(SIGSEGV, sigdie);
  signal(SIGXCPU, sigdie);
}

int main(int argc, char *argv[]) {
  try {
    setup_signal_handlers();
    parse_args(argc, argv);
  
    srandom(time(0));
  
    scanf("%d%d",&n,&m);
    R = matrix<double>(n, m+1);
    adj = matrix<double>(m+1, m+1);
  
    FOR(i,n) {
      FOR(j,m) {
        scanf("%lf",&R[i][j+1]);
      }
    }
  
    // product 0 is the null product
    ++m;
    FOR(i,n) R[i][0] = 0;
  
    FOR(i,n) {
      double cs;
      scanf("%lf",&cs);
  
      FOR(j,m) {
        R[i][j] = max(0.0, R[i][j] - cs);
      }
    }
  
    FOR(i,n) {
      scanf("%lf",&N[i]);
    }
  
    scanf("%lf",&delta); // ignore delta for now
    if (fabs(delta) > eps) {
      printf("WARNING: delta /= 0 and is ignored\n");
    }
  
    begin_time = get_cpu_usecs();
    int count = 0;

    do_heuristic_init();
    if (do_good) {
      do_good_heuristics(&count);
    } else {
      do_heuristics(&count);
    }

    assert((int)asst.size() == n);
    verify();
    if (do_dump_soln) dumpsoln();
    if (do_dump_pi) dumppi();
    if (do_dump_dir) dumpdir();
    last_words();
    return 0;
  } catch(std::bad_alloc &) {
    printf("Dying due to std::bad_alloc\n");
    last_words();
    fflush(stdout);
    return -1;
  }
}

namespace asdf {

struct dobson_soln_info {
  double profit;
  vector<int> parents;
  vector<double> prices;
  vector<int> buy;
};

void publish(const dobson_soln_info &soln) {
  FOR(i, m) pi[i] = soln.prices[i];
  FOR(i, n) asst[i] = soln.buy[i];
}

template <typename T>
struct queueue {
  vector<T> q;
  int q0, q1;
  queueue() : q0(0), q1(0) {}
  inline size_t size() const { return q1-q0; }
  inline void push(const T &t) { q.push_back(t); q1++; }
  inline void pop() { q0++; }
  inline const T &front() const { return q[q0]; }
  inline T &front() { return q[q0]; }
};

// When we run the Dobson-Kalish heuristic, there are m very closely related
// shortest-path problems to solve.  So closely related, in fact, that
// their adjacency matrices differ in exactly four rows and four columns.
// Rather than having to store t copies of the matrix (if we have t threads),
// we instead store one copy of the matrix, and special-case the rows and
// columns that are different.  (This is the sort of place where a reference-
// counting pointer would be incredibly useful, too.)
// Sadly, I never got around to actually implementing Dobson-Kalish in
// parallel, so this code may or may not work for that purpose.
// (Looking at it, I don't think it should.)
template <typename T> class dkmatrix {
 private:
  vector<vector<T>*> v;
  vector<int> badcols;
  int nbadcol;
  vector<vector<T> > thecols;
  int owner;
  int n;
  T dfl;
 public:
  dkmatrix(int n, T dfl) {
    v = vector<vector<T>*>(n);
    for (int i = 0; i < n; i++)
      v[i] = new vector<T>(n, dfl);
    nbadcol = 0;
    owner = 1;
    this->n = n;
    this->dfl = dfl;
  }
  dkmatrix(int n) {
    v = vector<vector<T>*>(n);
    for (int i = 0; i < n; i++)
      v[i] = new vector<T>(n, T());
    owner = 1;
    nbadcol = 0;
    this->n = n;
    dfl = T();
  }
  dkmatrix(const dkmatrix&foo) {
    v = foo.v;
    n = foo.n;
    nbadcol = 0;
    owner = 0;
  }
  ~dkmatrix() {
    if (owner) {
      for (int i = 0; i < n; i++)
        delete v[i];
    }
    else {
      for (int i = 0; i < nbadcol; i++)
        delete v[badcols[i]];
    }
  }
  private:
  inline __attribute__((always_inline)) vector<T> &operator[](int i) {
    for (int k = 0; k < nbadcol; k++)
      if (badcols[k] == i)
        return thecols[k];
    return *(v[i]);
  }
  inline __attribute__((always_inline))
  const vector<T> &operator[](int i) const {
    for (int k = 0; k < nbadcol; k++)
      if (badcols[k] == i)
        return thecols[k];
    return *(v[i]);
  }
  public:
  inline __attribute__((always_inline))
  T &at(int i, int j) { return (*this)[i][j]; }
  inline __attribute__((always_inline))
  const T &at(int i, int j) const { return (*this)[i][j]; }

  inline __attribute__((always_inline))
  void save_vertex(int i) {
    badcols.push_back(i);
    thecols.push_back(vector<T>(n, dfl));
    nbadcol++;
    v[i] = new vector<T>(n, dfl);
  }
};

class ShortestPaths {
  // Edge weights in the graph through which we are to find shortest paths.
  // We use this bizarre matrix class because I once had the idea to make
  // the Dobson-Kalish code use the other 63 processors in pilatus.
  dkmatrix<double> weights;
  // Shortest path tree.
  vector<int> parent;
  // How much do we care about the path from vertex i?
  vector<double> vertexwt;
  // How much does product i cost?
  vector<double> bestwt;
  int n;
 public:
  ShortestPaths(int n) : weights(n, 1e20) {
    this->n = n;
    parent = vector<int>(n, -1);
    vertexwt = vector<double>(n, 0);
    bestwt = vector<double>(n, 1e20);
    bestwt[0] = 0;
  }
  ShortestPaths(const ShortestPaths&sp) : weights(sp.weights) {
    n = sp.n;
    parent = sp.parent;
    vertexwt = sp.vertexwt;
    bestwt = sp.bestwt;
  }
  void save_vertex(int v) {
    weights.save_vertex(v);
  }
  // Why do these four methods exist?  They should invalidate some kind of
  // shortest-path length as well as setting the appropriate quantity, but
  // this winds up being incredibly fiddly --- we do this when we need to
  // instead.
  inline void set_edge(int a, int b, double w) {
    weights.at(a,b) = w;
  }
  inline double get_edge(int a, int b) {
    return weights.at(a,b);
  }
  inline void set_vertex_wt(int v, double w) {
    vertexwt[v] = w;
  }
  inline double get_vertex_wt(int v) {
    return vertexwt[v];
  }
  // Solve the shortest-path problem from scratch.
  void doit() {
    vector<int> dad(n,0);
    vector<int> hit(n,0);
    bestwt = vector<double>(n, 1e20);
    bestwt[0] = 0;
    dad[0] = 0;
    queueue<pair<double,int> > q;
    q.push(make_pair(0,0));
    while (q.size()) {
      int v = q.front().second; q.pop();
      for (int i = 0; i < n; i++)
        if (bestwt[i] > bestwt[v] + weights.at(i,v))
          bestwt[i] = bestwt[v] + weights.at(i,v),
          dad[i] = v,
          q.push(make_pair(bestwt[i],i));
    }
#ifdef DEBUG_SHORTEST_PATHS
    {
      vector<vector<double> > foo(n, vector<double>(n));
      for (int i=0;i<n;i++) {
        for (int j=0;j<n;j++)
          foo[i][j] = weights.at(i,j);
        foo[i][i] = 0;
      }
      for (int k=0;k<n;k++)
        for (int i=0;i<n;i++)
          for (int j=0;j<n;j++)
            foo[i][j] = min(foo[i][j], foo[i][k] + foo[k][j]);
      //for (int i=0;i<n;i++) printf("%lf\n", foo[i][i]);
  
      for (int i=0;i<n;i++)
        for (int j=0;j<n;j++)
          assert(bestwt[j] <= bestwt[i] + weights.at(j,i));

      FOR(i, n) if (i) {
        /*
        printf("%i %i:  %.0f %.0f %.0f %.0f\n", i, dad[i],
            bestwt[i], bestwt[dad[i]], weights.at(i, dad[i]),
            bestwt[dad[i]] + weights.at(i, dad[i]));
        */
        assert(bestwt[i] == bestwt[dad[i]] + weights.at(i, dad[i]));
      }
    }
#endif
    parent = dad;
  }
  // Update shortest path lengths for v and all descendants.
  // Assumes that, for each vertex w that is not v or a descendant of v,
  // bestwt[w] is larger than or equal to the length of some path from the
  // root to w.
  void redoit(int v) {
    queueue<int> q;
    int root = parent[v];
    if (v) bestwt[v] = 1e20;
    q.push(v);
    vector<int> kids;
    // Set path lengths for all kids of v to infinity.
    // Shortest paths to stuff other than kids of v will only go down, since
    // we've only screwed with edges involving v and parent[v].  However,
    // shortest paths to kids of v can go up, so we have to reset them.
    // This could probably be made faster by setting the distance to v and
    // each descendant to a better upper bound.
    while (q.size()) {
      int v = q.front(); q.pop();
      kids.push_back(v);
      for (int i = 1; i < n; i++) if (parent[i] == v)
        bestwt[i] = 1e20,
        q.push(i);
    }
    // So now bestwt[v] is an upper bound on the length of the shortest path
    // from n-1 to v for all vertices v.  bestwt[v] shall always have this
    // property; the purpose of the hacked Bellman-Ford-Moore below is to
    // make bestwt[v] equal to the length of the shortest path from n-1 to
    // v for all v.
    queueue<int> pq;
    vector<int> inqueue(n,0);
    pq.push(root);
    while (pq.size()) {
      int v = pq.front(); pq.pop();
      inqueue[v] = 0;
      // Yes, this is necessary.  Again, a small speedup could probably be
      // achieved by being careful about what bestwt[v] is before we run
      // Bellman-Ford-Moore.
      for (int i = 0; i < n; i++)
        if (bestwt[v] > bestwt[i] + weights.at(v,i))
          bestwt[v] = bestwt[i] + weights.at(v,i),
          parent[v] = i;
      // Now do paths to other vertices.  Make it so that we consider every
      // vertex to which we learnt of a shorter path.
      for (int i = 0; i < n; i++)
        if (bestwt[i] > bestwt[v] + weights.at(i,v)) {
          bestwt[i] = bestwt[v] + weights.at(i,v),
          parent[i] = v;
          if (!inqueue[i]) pq.push(i), inqueue[i] = 1;
        }
    }
#if DEBUG_SHORTEST_PATHS
    // Verify that the shortest path lengths are legitimate.
    {
      for (int i=0;i<n;i++)
        for (int j=0;j<n;j++)
          if (bestwt[j] > bestwt[i] + weights.at(j,i)) {
            printf("error %i %i\n", i, j);
            printf("kids were\n");
            for (int k=0;k<kids.size();k++)
              printf("%i ", kids[k]);
            printf("\n");
            exit(0);
          }
      FOR(i, n) if (i)
        assert(bestwt[i] == bestwt[parent[i]] + weights.at(i, parent[i]));
    }
#endif
  }
  double getprofit() {
    double p = 0;
    for (int i = 0; i < n; i++)
      p += vertexwt[i] * bestwt[i];
    return p;
  }
  int getparent(int v) { return parent[v]; }
  vector<int> getparents() { return parent; }
  double getprice(int v) { return bestwt[v]; }
  vector<double> getprices() { return bestwt; }
};

// Find the best prices for a given assignment of customers to products,
// and store all relevant information in the returned dobson_soln_info.
dobson_soln_info dobson_2(const vector<int>&buy) {
  ShortestPaths graph(m);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      graph.set_edge(i,j,1e20);
  for (int i = 0; i < n; i++)
    for (int k = 0; k < m; k++)
      graph.set_edge(buy[i],k,min(graph.get_edge(buy[i],k),
                                  R[i][buy[i]]-R[i][k]));
  for (int i = 0; i < n; i++)
    graph.set_vertex_wt(buy[i], graph.get_vertex_wt(buy[i]) + N[i]);
  graph.doit();
  dobson_soln_info rv;
  rv.profit = graph.getprofit();
  rv.parents = graph.getparents();
  rv.prices = graph.getprices();
  return rv;
}

// Find the Dobson-Kalish move that gives the greatest increase in the
// objective function.  This function does most of the real work in
// the Dobson-Kalish heuristic.
int dobson_getbestmove(
        vector<int> buy, const vector<int>&parents, double oldprofit,
        double &nuprofit) {
  ShortestPaths graph(m);
  // set each edge to the appropriate value.
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      graph.set_edge(i,j,1e20);
  for (int i = 0; i < n; i++)
    for (int k = 0; k < m; k++)
      graph.set_edge(buy[i],k,min(graph.get_edge(buy[i],k),
                                  R[i][buy[i]]-R[i][k]));
  for (int i = 0; i < n; i++)
    graph.set_vertex_wt(buy[i], graph.get_vertex_wt(buy[i]) + N[i]);

  // solve the shortest-path problem.
  graph.doit();

  int bestmove = -1;
  double bestprofit = oldprofit;
  // bestprice[j] is set to the highest price that a product can command
  // under the current assignment customers.
  vector<double> bestprice(m, 1e20);
  for (int i = 0; i < n; i++) //if (buy[i] != 0)
    bestprice[buy[i]] = min(bestprice[buy[i]],
                            R[i][buy[i]] - R[i][parents[buy[i]]]);

  // For each customer i who constrains the increase of the price of some
  // product...
  for (int i = 0; i < n; i++)
   if (   buy[i] != 0
       && R[i][buy[i]] - R[i][parents[buy[i]]] == bestprice[buy[i]]) {
    // Assign i to the parent product.
    int old = buy[i];
    buy[i] = parents[buy[i]];

    // Construct the new graph.  (We do this in an incredibly brutal way;
    // the old graph is not actually saved, only the shortest path lengths.)
    // This should really be fixed so that only the edges that change are
    // updated.
    ShortestPaths graph2(graph);
    graph2.save_vertex(old);
    graph2.save_vertex(buy[i]);
    for (int j = 0; j < m; j++)
      graph2.set_edge(old,j,1e20),
      graph2.set_edge(j,old,1e20),
      graph2.set_edge(buy[i],j,1e20),
      graph2.set_edge(j,buy[i],1e20);

    for (int j = 0; j < n; j++)
      if (buy[j] == old || (buy[j] == buy[i] && buy[j] != 0)) {
        for (int k = 0; k < m; k++)
          graph2.set_edge(buy[j],k,min(graph2.get_edge(buy[j],k),
                                       R[j][buy[j]] - R[j][k]));
      }
      else if (buy[j] != 0) {
        graph2.set_edge(buy[j], buy[i],
                        min(graph2.get_edge(buy[j], buy[i]),
                            R[j][buy[j]] - R[j][buy[i]]));
        graph2.set_edge(buy[j], old,
                        min(graph2.get_edge(buy[j], old),
                            R[j][buy[j]] - R[j][old]));
      }
    graph2.set_vertex_wt(old, graph2.get_vertex_wt(old) - N[i]);
    graph2.set_vertex_wt(buy[i], graph2.get_vertex_wt(buy[i]) + N[i]);
    // Solve the new shortest path problem using the information from the
    // old shortest path problem, and figure out whether we did the right
    // thing.
    graph2.redoit(old);
    double prof = graph2.getprofit();
    graph2.set_vertex_wt(old, graph2.get_vertex_wt(old) + N[i]);
    graph2.set_vertex_wt(buy[i], graph2.get_vertex_wt(buy[i]) - N[i]);

    if (prof > bestprofit)
      bestprofit = prof, bestmove = i;
    buy[i] = old;
  }
  nuprofit = bestprofit;
  return bestmove;
}

}

// Do ONE iteration of DK.
bool torsdk() {
  vector<int> old_asst = asst;

  asdf::dobson_soln_info foo = asdf::dobson_2(asst);
  double oldprofit = foo.profit, nuprofit;
  int bestmove =
      asdf::dobson_getbestmove(asst, foo.parents, oldprofit, nuprofit);
  pi = foo.prices;
  if (bestmove == -1) {
    printf("Failed to do a DK move.\n");
    return false;
  }

  printf("  Bumped %i from %i to %i changing profit from %f to %f\n",
          bestmove, asst[bestmove], foo.parents[asst[bestmove]],
          oldprofit, nuprofit);
  asst[bestmove] = foo.parents[asst[bestmove]];
  return true;
}
