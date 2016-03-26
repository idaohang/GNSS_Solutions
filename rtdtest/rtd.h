#ifndef _RTD_H_
#define _RTD_H_

#include "src\rtklib.h"

#define SOL_SUCCESS 1
#define SOL_FAIL 0

#define WGS84_a       6.3781370e6  //长半轴
#define WGS84_GM      3.986004418e14  //地球引力常数（含大气层）
#define WGS84_omegae  7.2921151467e-5  //地球自转角速度
#define WGS84_ell     298.257223563  //扁率
#define WGS84_U0      62636860.8497  //椭球正常重力位
#define WGS84_gammae  9.7803267714  //赤道正常重力位
#define WGS84_b       6356752.3142  //短半轴
#define WGS84_C20     -484.16685e-6  //引力位二阶谐系数
#define WGS84_e2      0.00669437999013  //第一偏心率平方
#define WGS84_e       sqrt(WGS84_e2)  //第一偏心率
#define WGS84_ee2     0.006739496742227  //第二偏心率平方
#define WGS84_ee      sqrt(WGS84_ee2)  //第二偏心率平方
#define F             -4.442807633e-10

typedef struct {
	int sys;
	int sat;
	int prn;
	double pos[3];      /* satellite position (m) (ecef) */
	double vel[3];      /* satellite velocity (m/s) (ecef) */
	double acc[3];      /* satellite acceleration (m/s^2) (ecef) */
	double af0, af1;     /* satellite clock-offset/drift (s,s/s) */
} sat_t;

#ifdef _cplusplus
extern "C" {
#endif
	//declare functions for rtd
	int open_rt_stream(char* ntrip_str);
	int run_rtd();
	int my_rtd(obs_t obss, nav_t navs, sol_t sol);

	//Satellite position velocity clock error clock rate

	void writeFile();

#ifdef _cplusplus
}
#endif

#endif