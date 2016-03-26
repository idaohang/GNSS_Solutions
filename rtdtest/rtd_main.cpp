#include <iostream>
#include <string>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iomanip>

#include "src/rtklib.h"
#include "rtd.h"

using namespace std;
FILE* g_mylog = NULL;
FILE* binFile = NULL;

#define MAXSTR      5                  /* max number of streams */
static obs_t obss = { 0 };          /* observation data */
static nav_t navs = { 0 };          /* navigation data */
int main(int argc, char* argv[])
{
	binFile = fopen("log/binfile.bin", "wb");
	g_mylog = fopen("mytrace.log", "w+");

	open_rt_stream(NULL);
	run_rtd();

	return 0;
}

static strsvr_t  mysvr_obs, mysvr_eph, mysvr_ssr;
static strconv_t *myconv_obs[MAXSTR] = { NULL }, *myconv_eph[MAXSTR] = { NULL }, *myconv_ssr[MAXSTR] = { NULL };
int open_rt_stream(char* ntrip_str) {

	FILE* ephfile = NULL, *ssrfile = NULL, *obsfile = NULL;
	FILE* recordfile = NULL;

	char *paths_obs[MAXSTR] = { NULL };
	char *paths_eph[MAXSTR] = { NULL };
	char *paths_ssr[MAXSTR] = { NULL };

	char* ntrip_obs = "001@whu:123|119.97.244.11:2102/VRS2";  //obs
	//char* ntrip_obs = "gxwang:wang123|58.49.58.149:2101/ABMF0";  //obs
	char* ntrip_eph = "gxwang:wang123|58.49.58.149:2101/AREG7";	//IGS  eph
	char* ntrip_ssr = "gxwang:wang123|58.49.58.149:2101/CLK91";	//IGS  ssr
	paths_obs[0] = ntrip_obs;  paths_obs[1] = "mytest.dat";
	paths_eph[0] = ntrip_eph;  paths_eph[1] = ":2100";
	paths_ssr[0] = ntrip_ssr;  paths_ssr[1] = "mytest.dat";	
	
	char satstr[3] = { 0 };
	int sys = -1, prn = -1;
	char timestr[20] = { 0 };
	obsd_t obs[MAXOBS];
	//double rb[3] = { 0.0 }, pos[3] = { 0.0 }, e[3] = { 0.0 };
	/*strsvr_t  mysvr_obs, mysvr_eph, mysvr_ssr;*/
	int  dispint = 1000;
	int opts[] = { 10000, 10000, 2000, 32768, 10, 0, 30 };
	double stapos[3] = { 0 };
	/*strconv_t *myconv_obs[MAXSTR] = { NULL }, *myconv_eph[MAXSTR] = { NULL }, *myconv_ssr[MAXSTR] = { NULL };*/
	int types[MAXSTR] = { STR_NTRIPCLI, STR_FILE };
	int types0[MAXSTR] = { STR_NTRIPCLI, STR_TCPSVR };
	
	char* msgout_obs = "1004(1),1012(1),1006(10),1008(10)";
	char* msgout_eph = "1007(10),1008(10),1019(15),1033(10),1045(15),1046(15),1077(1),1087(1),1097(1),1107(1),1127(1)";
	char* msgout_ssr = "1059(5),1060(5),1065(5),1066(5)";
	char* opt = "";
	int test_svr[3] = { 0 };
	int strstat[3][MAXSTR] = { { 0 } };
	int strbyte[3][MAXSTR] = { { 0 } };
	int bps[3][MAXSTR] = { { 0 } };
	char strmsg[3][MAXSTRMSG] = { {""} };
	const char ss[] = { 'E', '-', 'W', 'C', 'C' };
	char buff[256], *p;
	int i = 0, j = 0, n = 1, t1;  // n的值至少是1
	int nepoch;

	obsfile = fopen("log/obs.log", "w+");
	ephfile = fopen("log/eph.log", "w+");
	ssrfile = fopen("log/ssr.log", "w+");
	recordfile = fopen("log/record.log", "w+");
	fflush(obsfile);  //全缓冲
	fflush(ephfile);
	fflush(ssrfile);
	fflush(recordfile);

	myconv_obs[0] = strconvnew(STRFMT_BINEX, STRFMT_RTCM3, msgout_obs, 0, 0, opt);
	myconv_eph[0] = strconvnew(STRFMT_RTCM3, STRFMT_RTCM3, msgout_eph, 0, 0, opt);
	myconv_eph[1] = strconvnew(STRFMT_RTCM3, STRFMT_RTCM3, msgout_eph, 0, 0, opt);
	myconv_ssr[0] = strconvnew(STRFMT_RTCM3, STRFMT_RTCM3, msgout_ssr, 0, 0, opt);

	strsvrinit(&mysvr_obs, n);
	strsvrinit(&mysvr_eph, n);
	strsvrinit(&mysvr_ssr, n);
	//启动服务线程
	test_svr[0] = strsvrstart(&mysvr_obs, opts, types, paths_obs, myconv_obs, NULL, stapos);
	if (!test_svr[0]) {
		printf("stream(RTCM) server start error\n");
		system("pause");
		return -1;
	}

	test_svr[1] = strsvrstart(&mysvr_eph, opts, types0, paths_eph, myconv_eph, NULL, stapos);
	if (!test_svr[1]) {
		printf("stream(RTCM) server start error\n");
		system("pause");
		return -1;
	}

	test_svr[2] = strsvrstart(&mysvr_ssr, opts, types, paths_ssr, myconv_ssr, NULL, stapos);
	if (!test_svr[2]) {
		printf("stream(RTCM) server start error\n");
		system("pause");
		return -1;
	}
}


int run_rtd() {
	int nepoch = 0, obsflag = 0, nobs = 0;
	int sys;
	char satstr[10];
	gtime_t time;
	int  dispint = 1000;
	char timestr[20];
	char strmsg[3][MAXSTRMSG] = { { "" } };
	unsigned int t1;
	int ephopt = EPHOPT_SSRAPC;
	double* rs = (double*)malloc(MAXOBS * 6);
	double* dts = (double*)malloc(MAXOBS * 2);
	double* var = (double*)malloc(MAXOBS);
	int* svh = (int*)malloc(MAXOBS);
	FILE * recordfile = fopen("log/recordfile.log", "w");
	fflush(recordfile);
	//开始循环计算
	while (true) {
		//锁定观测值接受线程
		lock(&mysvr_obs.lock);
		obsflag = myconv_obs[0]->out.obsflag;
		if (obsflag == 1) {  // obs data complete flag 1:ok
			obss = myconv_obs[0]->out.obs;  // update static variable
			navs = myconv_obs[0]->out.nav;
			for (int i = 0; i < obss.n; i++) {
				obss.data[i].rcv = 1;  // 1 表示基准站
			}
			
			nepoch = sortobs(&obss); //sort and unique observation data by time, rcv, sat
			time = obss.data->time;
			time2str(time, timestr, 3);
		}

		unlock(&mysvr_obs.lock);

		//锁定星历接受线程
		lock(&mysvr_eph.lock);
		//navs = myconv_eph[0]->out.nav;  //update the static variable
		unlock(&mysvr_eph.lock);

		//锁定ssr改正数线程
		lock(&mysvr_ssr.lock);
		memcpy(navs.ssr, myconv_ssr[0]->rtcm.ssr, sizeof(ssr_t)*MAXSAT);
		unlock(&mysvr_ssr.lock);

		//获得了一个历元解算需要的数据，进行解算
		t1 = tickget();
		if (obsflag==1)
		{
			cout << "one epoch" << endl;
		}
		if (obsflag == 1&&false) {
			cout << "one epoch" << endl;
			nobs = obss.n;
			printf("%s  %d\n", timestr, nepoch++);

			fprintf(recordfile, "%s\n", timestr);
			satposs(time, obss.data, obss.n, &navs, EPHOPT_BRDC, rs, dts, var, svh);
			fprintf(recordfile, "%s\n", "Without ssr\n");
			for (int i = 0; i < nobs; i++) {
				sys = satsys(obss.data[i].sat, NULL);
				if (sys != SYS_GPS) {
					continue;
				}
				satno2id(obss.data[i].sat, satstr);
				fprintf(recordfile, "%s: X:%lf; Y:%lf; Z:%lf\n", satstr, rs[6 * i + 0], rs[6 * i + 1], rs[6 * i + 2]);

				fwrite(&(obss), sizeof(obs_t), 3, binFile);
				//_int64 time64 = obss.data[i].time.time;
				//fwrite(&time64, sizeof(_int64), 1, binFile);

				//int num = 65450;
				//fwrite(&num, sizeof(int), 1, binFile);
			}

			satposs(time, obss.data, obss.n, &navs, EPHOPT_SSRAPC, rs, dts, var, svh);
			fprintf(recordfile, "%s\n", "With ssr\n");
			for (int i = 0; i < nobs; i++) {
				sys = satsys(obss.data[i].sat, NULL);
				if (sys != SYS_GPS) {
					continue;
				}
				satno2id(obss.data[i].sat, satstr);
				fprintf(recordfile, "%s: X:%lf; Y:%lf; Z:%lf\n", satstr, rs[6 * i + 0], rs[6 * i + 1], rs[6 * i + 2]);
			}
		}
		
		sleepms(dispint - (tickget() - t1));
	}
}