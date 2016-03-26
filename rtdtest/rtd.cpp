#include <iostream>
#include "src\rtklib.h"
#include "rtd.h"
using namespace std;

int traverseNavs(const obs_t* obss, const nav_t* navs);
int my_rtd(obs_t* obss, nav_t* navs, sol_t sol) {
	if (traverseNavs(obss, navs) == 0) {
		return SOL_FAIL;
	}
	return SOL_FAIL;
}

/*判断观测数据中的所有卫星健康度
  return 1：sat count > 6
  return 0: sat count <= 6
*/
int traverseNavs(const obs_t* obss, const nav_t* navs) {
	int satCount = obss->n;
	for (int i = 0; i < obss->n; i++) {
		int sat = obss->data[i].sat;
		int svh = navs->eph[sat - 1].svh;
		if (satexclude(sat, svh, NULL) == 1) {
			satCount--;
			continue;
		}
	}
	if (satCount < 6) {
		return 0;
	} else {
		return 1;
	}
}
