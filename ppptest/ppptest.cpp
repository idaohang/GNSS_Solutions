// ppptest.cpp : 定义控制台应用程序的入口点。
//

#include "src/rtklib.h"

#include <iostream>
#include <string>
using namespace std;
#include <math.h>
#include <fstream>
#include <iomanip>
//
////全局变量
//extern pcvs_t pcvss;        /* receiver antenna parameters */
//extern pcvs_t pcvsr;        /* satellite antenna parameters */
//extern obs_t obss;          /* observation data */
//extern nav_t navs;          /* navigation data */
//extern sbs_t sbss;          /* sbas messages */
//extern lex_t lexs;          /* lex messages */
//extern sta_t stas[MAXRCV];      /* station infomation */
//extern int nepoch;            /* number of observation epochs */
//extern int iobsu ;            /* current rover observation data index */
//extern int iobsr ;            /* current reference observation data index */
//extern int isbs  ;            /* current sbas message index */
//extern int ilex  ;            /* current lex message index */
//extern int revs  ;            /* analysis direction (0:forward,1:backward) */
//extern int aborts;            /* abort status */
//extern sol_t *solf;             /* forward solutions */
//extern sol_t *solb;             /* backward solutions */
//extern double *rbf;             /* forward base positions */
//extern double *rbb;             /* backward base positions */
//extern int isolf;             /* current forward solutions index */
//extern int isolb;             /* current backward solutions index */
//extern char proc_rov [64];   /* rover for current processing */
//extern char proc_base[64];   /* base station for current processing */
//extern char rtcm_file[1024]; /* rtcm data file */
//extern char rtcm_path[1024]; /* rtcm data path */
//extern rtcm_t rtcm;             /* rtcm control struct */
//extern FILE *fp_rtcm;      /* rtcm data file pointer */

FILE* g_mylog = NULL;
static volatile int intrflg = 0;         /* interrupt flag */
#define MAXSTR      5                  /* max number of streams */

//#define  MAXBUFLEN 5120 
//#define PORT 8888
//#define DEST_IP_ADDR "127.0.0.1" //想要链接的目标Server address  

//
//typedef struct
//{        /* stream server type */
//	int state;          /* server state (0:stop,1:running) */
//	int cycle;          /* server cycle (ms) */
//	int buffsize;       /* input/monitor buffer size (bytes) */
//	int nmeacycle;      /* NMEA request cycle (ms) (0:no) */
//	int nstr;           /* number of streams (1 input + (nstr-1) outputs */
//	int npb;            /* data length in peek buffer (bytes) */
//	double nmeapos[3];  /* NMEA request position (ecef) (m) */
//	unsigned char *buff; /* input buffers */
//	unsigned char *pbuf; /* peek buffer */
//	unsigned int tick;  /* start tick */
//	stream_t stream[16]; /* input/output streams */
//	strconv_t *conv[16]; /* stream converter */
//	thread_t thread;    /* server thread */
//	lock_t lock;        /* lock flag */
//} mySVR_t;
//




/*解码和芯星通数据*/
int getUB380Data(char* ntrip_str, int dispint)
{

	if (ntrip_str == NULL)  { ntrip_str = "001@whu:123|119.97.244.11:2101/WHU01"; }  //"gxwang:wang123@58.49.58.149:2101/CLK91"
	if (dispint == 0)   { dispint = 1000; }  // 5 seconds by default
	//for the real time test
	strsvr_t  mysvr;
	int opts[] = { 10000, 10000, 2000, 32768, 10, 0, 30 };
	double stapos[3] = { 0 };
	//char*  myaddress = "yangmf:123456@59.175.223.165:2101/WHU01";
	char *paths[MAXSTR] = { "001@whu:123|119.97.244.11:2101/ZY01" };
	paths[0] = ntrip_str;
	//paths[0] = "gxwang:wang123@58.49.58.149:2101/CLK20";
	paths[1] = "mytest.dat";
	int types[MAXSTR] = { STR_NTRIPCLI, STR_FILE };
	strconv_t *myconv[MAXSTR] = { NULL };
	char *msgout = "1060,1066,1238";
	char* opt = "";
	int strstat[MAXSTR] = { 0 };
	int strbyte[MAXSTR] = { 0 };
	int bps[MAXSTR] = { 0 };
	char strmsg[MAXSTRMSG] = "";
	const char ss[] = { 'E', '-', 'W', 'C', 'C' };
	char buff[256], *p;
	int i = 0, j = 0, n = 1;  // n的值至少是1
	myconv[0] = strconvnew(STRFMT_UB380, STRFMT_RTCM3, msgout, 0, 0, opt);
	fstream  recodfile("log/igs.rec", ios::out);
	nav_t mynav = { { 0 } };
	obs_t myobs = { { 0 } };
	strsvrinit(&mysvr, n);

	int test_svr = strsvrstart(&mysvr, opts, types, paths, myconv, NULL, stapos);
	if (!test_svr)
	{
		printf("stream server start error\n");
		return -1;
	}

	gtime_t  mytime[32] = { { 0 } };
	char timestr[20] = { 0 }, satstr[3] = { 0 };
	int tag = 0, prn = -1, sys = -1;

	rnxopt_t rnxoption;
	rnxoption.navsys = SYS_GPS;
	rnxoption.rnxver = 3.0;
	rnxoption.nobs[0] = rnxoption.nobs[1] = rnxoption.nobs[2] = rnxoption.nobs[3] = rnxoption.nobs[4] = rnxoption.nobs[5] = MAXOBS;
	rnxoption.obstype = OBSTYPE_PR | OBSTYPE_CP | OBSTYPE_SNR | OBSTYPE_DOP;
	rnxoption.freqtype = FREQTYPE_L1 | FREQTYPE_L2 | FREQTYPE_L5 | FREQTYPE_L7;

	FILE* rnxfile = fopen("log/rnxobs.obs", "w+");
	if (rnxfile == NULL)
	{
		printf("记录文件打开错误!\n");
		exit(0);
	}

	for (intrflg = 0; !intrflg;)
	{

		lock(&mysvr.lock);

		//printf("%d\n",myconv[0]->raw.ephsat);
		//memcpy(&mynav,&myconv[0]->raw.nav,sizeof(nav_t));
		mynav = myconv[0]->raw.nav;
		myobs = myconv[0]->raw.obs;
		printf("星历信息:%d---%d\n", myconv[0]->raw.ephsat, mynav.n);
		for (i = 0; i < mynav.n; i++)
		{
			sys = satsys(mynav.eph[i].sat, &prn);
			if (sys != SYS_GPS)
			{
				continue;
			}
			satno2id(mynav.eph[i].sat, satstr);
			printf("%s:%.4f ", satstr, mynav.eph[i].A);
		}
		printf("\n");


		time2str(myobs.data[0].time, timestr, 3);
		printf("观测值信息:\n");
		printf("epochtime:%s             \n", timestr);
		for (i = 0; i < myobs.n; i++)
		{
			sys = satsys(myobs.data[i].sat, &prn);
			if (sys != SYS_GPS)
			{
				continue;
			}
			satno2id(myobs.data[i].sat, satstr);
			printf("sat: %s obs:%.3f  %.3f %.3f %.3f \n", satstr,
				myobs.data[i].P[0], myobs.data[i].P[1],
				myobs.data[i].L[0], myobs.data[i].L[1]);
		}
		printf("\n");

		for (i = 0; i < myobs.n; i = j)
		{
			for (j = i; j<myobs.n; j++)
			{
				if (timediff(myobs.data[j].time, myobs.data[i].time)>0.001) break;
			}
			outrnxobsb(rnxfile, &rnxoption, myobs.data + i, j - i, 0);
		}

		//if( tag == 0 )
		//{
		//	for( i =0 ; i<32; i++ )
		//	{
		//		mytime[i] = myconv[0]->rtcm.ssr[i].t0[0];
		//	}
		//	tag = 1;
		//}
		//else
		//{
		//	if( fabs(timediff( myconv[0]->rtcm.ssr[i].t0[0],mytime[i])) > 0.01   )
		//	{
		//		recodfile<<timestr<<" ";
		//	}
		//	for( i = 0 ; i< 32; i++ )
		//	{
		//		if( fabs(timediff( myconv[0]->rtcm.ssr[i].t0[0],mytime[i])) > 0.01   )
		//		{
		//			time2str(myconv[0]->rtcm.ssr[i].t0[0],timestr,3);
		//			printf("星历参考时刻:%s\n",timestr);
		//			printf("卫星号：%d --> Along:%4.3f Radial:%4.3f Cross:%4.3f Along_dot:%8.7f Radial_dot:%8.7f Cross_dot:%8.7f\n",
		//				i+1,myconv[0]->rtcm.ssr[i].deph[0],myconv[0]->rtcm.ssr[i].deph[1],myconv[0]->rtcm.ssr[i].deph[2],
		//				myconv[0]->rtcm.ssr[i].ddeph[0],myconv[0]->rtcm.ssr[i].ddeph[1],myconv[0]->rtcm.ssr[i].ddeph[2]);
		//			printf("钟参数:-->C0:%4.3f C1:%.4e C2:%.4e\n",myconv[0]->rtcm.ssr[i].dclk[0],myconv[0]->rtcm.ssr[i].dclk[1],myconv[0]->rtcm.ssr[i].dclk[2]);

		//			recodfile<<setw(8)<< myconv[0]->rtcm.ssr[i].deph[0]<<" "; //<<myconv[0]->rtcm.ssr[i].deph[1]<<" "myconv[0]->rtcm.ssr[i].deph[2]<< 

		//			mytime[i] =  myconv[0]->rtcm.ssr[i].t0[0];
		//			if( i == 31)  recodfile<<endl; 	
		//		}
		//	}
		//}




		/*double dclk = myconv[0]->rtcm.ssr[0].dclk[0];
		printf("*****%.4f*****    \n",dclk);*/
		unlock(&mysvr.lock);
		/* get stream server status */
		strsvrstat(&mysvr, strstat, strbyte, bps, strmsg);

		/* show stream server status */
		for (i = 0, p = buff; i < 1; i++)
		{
			p += sprintf(p, "%c", ss[strstat[i] + 1]);
		}

		printf("%s [%s] %10d B %7d bps %s\n",
			time_str(utc2gpst(timeget()), 0), buff, strbyte[0], bps[0], strmsg);

		sleepms(dispint);
	}


}

/*binex数据测试*/
int getBinexData(char* ntrip_str, int dispint)
{
	if (ntrip_str == NULL)  { ntrip_str = "001@whu:123|119.97.244.11:2101/GK03"; }  //"gxwang:wang123@58.49.58.149:2101/CLK91"
	if (dispint == 0)   { dispint = 1000; }  // 5 seconds by default
	//for the real time test
	strsvr_t  mysvr;
	int opts[] = { 10000, 10000, 2000, 32768, 10, 0, 30 };
	double stapos[3] = { 0 };
	//char*  myaddress = "yangmf:123456@59.175.223.165:2101/WHU01";
	char *paths[MAXSTR] = { "001@whu:123|119.97.244.11:2101/ZY01" };
	paths[0] = ntrip_str;
	//paths[0] = "gxwang:wang123@58.49.58.149:2101/CLK20";
	paths[1] = "mytest.dat";
	int types[MAXSTR] = { STR_NTRIPCLI, STR_FILE };
	strconv_t *myconv[MAXSTR] = { NULL };
	char *msgout = "1060,1066,1238";
	char* opt = "";
	int strstat[MAXSTR] = { 0 };
	int strbyte[MAXSTR] = { 0 };
	int bps[MAXSTR] = { 0 };
	char strmsg[MAXSTRMSG] = "";
	const char ss[] = { 'E', '-', 'W', 'C', 'C' };
	char buff[256], *p;
	int i = 0, j = 0, n = 1;  // n的值至少是1
	myconv[0] = strconvnew(STRFMT_BINEX, STRFMT_RTCM3, msgout, 0, 0, opt);
	fstream  recodfile("log/igs.rec", ios::out);
	nav_t mynav = { { 0 } };
	obs_t myobs = { { 0 } };
	strsvrinit(&mysvr, n);

	int test_svr = strsvrstart(&mysvr, opts, types, paths, myconv, NULL, stapos);
	if (!test_svr)
	{
		printf("stream server start error\n");
		return -1;
	}

	gtime_t  mytime[32] = { { 0 } };
	char timestr[20] = { 0 }, satstr[3] = { 0 };
	int tag = 0, prn = -1, sys = -1;
	for (intrflg = 0; !intrflg;)
	{
		lock(&mysvr.lock);
		//printf("%d\n",myconv[0]->raw.ephsat);
		//memcpy(&mynav,&myconv[0]->raw.nav,sizeof(nav_t));
		mynav = myconv[0]->raw.nav;
		myobs = myconv[0]->raw.obs;
		printf("星历信息:%d---%d\n", myconv[0]->raw.ephsat, mynav.n);
		for (i = 0; i < mynav.n; i++)
		{
			sys = satsys(mynav.eph[i].sat, &prn);
			if (sys != SYS_GPS)
			{
				continue;
			}
			satno2id(mynav.eph[i].sat, satstr);
			printf("%s:%.4f ", satstr, mynav.eph[i].A);
		}
		printf("\n");

		time2str(myobs.data[0].time, timestr, 3);
		printf("观测值信息:\n");
		printf("epochtime:%s             \n", timestr);
		for (i = 0; i < myobs.n; i++)
		{
			sys = satsys(myobs.data[i].sat, &prn);
			if (sys != SYS_GPS)
			{
				continue;
			}
			satno2id(myobs.data[i].sat, satstr);
			printf("sat: %s obs:%.3f  %.3f %.3f %.3f \n", satstr,
				myobs.data[i].P[0], myobs.data[i].P[1],
				myobs.data[i].L[0], myobs.data[i].L[1]);
		}
		printf("\n");

		//if( tag == 0 )
		//{
		//	for( i =0 ; i<32; i++ )
		//	{
		//		mytime[i] = myconv[0]->rtcm.ssr[i].t0[0];
		//	}
		//	tag = 1;
		//}
		//else
		//{
		//	if( fabs(timediff( myconv[0]->rtcm.ssr[i].t0[0],mytime[i])) > 0.01   )
		//	{
		//		recodfile<<timestr<<" ";
		//	}
		//	for( i = 0 ; i< 32; i++ )
		//	{
		//		if( fabs(timediff( myconv[0]->rtcm.ssr[i].t0[0],mytime[i])) > 0.01   )
		//		{
		//			time2str(myconv[0]->rtcm.ssr[i].t0[0],timestr,3);
		//			printf("星历参考时刻:%s\n",timestr);
		//			printf("卫星号：%d --> Along:%4.3f Radial:%4.3f Cross:%4.3f Along_dot:%8.7f Radial_dot:%8.7f Cross_dot:%8.7f\n",
		//				i+1,myconv[0]->rtcm.ssr[i].deph[0],myconv[0]->rtcm.ssr[i].deph[1],myconv[0]->rtcm.ssr[i].deph[2],
		//				myconv[0]->rtcm.ssr[i].ddeph[0],myconv[0]->rtcm.ssr[i].ddeph[1],myconv[0]->rtcm.ssr[i].ddeph[2]);
		//			printf("钟参数:-->C0:%4.3f C1:%.4e C2:%.4e\n",myconv[0]->rtcm.ssr[i].dclk[0],myconv[0]->rtcm.ssr[i].dclk[1],myconv[0]->rtcm.ssr[i].dclk[2]);

		//			recodfile<<setw(8)<< myconv[0]->rtcm.ssr[i].deph[0]<<" "; //<<myconv[0]->rtcm.ssr[i].deph[1]<<" "myconv[0]->rtcm.ssr[i].deph[2]<< 

		//			mytime[i] =  myconv[0]->rtcm.ssr[i].t0[0];
		//			if( i == 31)  recodfile<<endl; 	
		//		}
		//	}
		//}

		/*double dclk = myconv[0]->rtcm.ssr[0].dclk[0];
		printf("*****%.4f*****    \n",dclk);*/
		unlock(&mysvr.lock);
		/* get stream server status */
		strsvrstat(&mysvr, strstat, strbyte, bps, strmsg);

		/* show stream server status */
		for (i = 0, p = buff; i < 1; i++)
		{
			p += sprintf(p, "%c", ss[strstat[i] + 1]);
		}

		printf("%s [%s] %10d B %7d bps %s\n",
			time_str(utc2gpst(timeget()), 0), buff, strbyte[0], bps[0], strmsg);

		sleepms(dispint);
	}
}

/*get the orbit and clock correction from ntrip rtcm(igs)*/
int getBrdcCorrection(char* ntrip_str, int dispint)
{
	if (ntrip_str == NULL)  { ntrip_str = "gxwang:wang123|58.49.58.149:2101/CLK91"; }  //"gxwang:wang123@58.49.58.149:2101/CLK91"
	if (dispint == 0)   { dispint = 1000; }  // 5 seconds by default
	//for the real time test
	strsvr_t  mysvr;
	int opts[] = { 10000, 10000, 2000, 32768, 10, 0, 30 };
	double stapos[3] = { 0 };
	//char*  myaddress = "yangmf:123456@59.175.223.165:2101/WHU01";
	char *paths[MAXSTR] = { "yangmf:123456@59.175.223.165:2101/WHU01" };
	paths[0] = ntrip_str;
	//paths[0] = "gxwang:wang123@58.49.58.149:2101/CLK20";
	paths[1] = "mytest.dat";
	int types[MAXSTR] = { STR_NTRIPCLI, STR_FILE };
	strconv_t *myconv[MAXSTR] = { NULL };
	char *msgout = "1060,1066,1238";
	char* opt = "";
	int strstat[MAXSTR] = { 0 };
	int strbyte[MAXSTR] = { 0 };
	int bps[MAXSTR] = { 0 };
	char strmsg[MAXSTRMSG] = "";
	const char ss[] = { 'E', '-', 'W', 'C', 'C' };
	char buff[256], *p;
	int i = 0, j = 0, n = 1;  // n的值至少是1
	myconv[0] = strconvnew(STRFMT_RTCM3, STRFMT_RTCM3, msgout, 0, 0, opt);
	fstream  recodfile("log/igs.rec", ios::out);
	strsvrinit(&mysvr, n);

	int test_svr = strsvrstart(&mysvr, opts, types, paths, myconv, NULL, stapos);
	if (!test_svr)
	{
		printf("stream server start error\n");
		return -1;
	}

	gtime_t  mytime[32] = { { 0 } };
	char timestr[20] = { 0 };
	int tag = 0;
	for (intrflg = 0; !intrflg;)
	{
		lock(&mysvr.lock);

		//if( /*myconv[0]->rtcm.ssr[i].dclk[0] != 0.0 &&*/ myconv[0]->rtcm.ssr[0].update == 1 )
		{
			if (tag == 0)
			{
				for (i = 0; i<32; i++)
				{
					mytime[i] = myconv[0]->rtcm.ssr[i].t0[0];
				}
				tag = 1;
			}
			else
			{
				if (fabs(timediff(myconv[0]->rtcm.ssr[i].t0[0], mytime[i])) > 0.01)
				{
					recodfile << timestr << " ";
				}
				for (i = 0; i< 32; i++)
				{
					if (fabs(timediff(myconv[0]->rtcm.ssr[i].t0[0], mytime[i])) > 0.01)
					{
						time2str(myconv[0]->rtcm.ssr[i].t0[0], timestr, 3);
						printf("星历参考时刻:%s\n", timestr);
						printf("卫星号：%d --> Along:%4.3f Radial:%4.3f Cross:%4.3f Along_dot:%8.7f Radial_dot:%8.7f Cross_dot:%8.7f\n",
							i + 1, myconv[0]->rtcm.ssr[i].deph[0], myconv[0]->rtcm.ssr[i].deph[1], myconv[0]->rtcm.ssr[i].deph[2],
							myconv[0]->rtcm.ssr[i].ddeph[0], myconv[0]->rtcm.ssr[i].ddeph[1], myconv[0]->rtcm.ssr[i].ddeph[2]);
						printf("钟参数:-->C0:%4.3f C1:%.4e C2:%.4e\n", myconv[0]->rtcm.ssr[i].dclk[0], myconv[0]->rtcm.ssr[i].dclk[1], myconv[0]->rtcm.ssr[i].dclk[2]);

						recodfile << setw(8) << myconv[0]->rtcm.ssr[i].deph[0] << " "; //<<myconv[0]->rtcm.ssr[i].deph[1]<<" "myconv[0]->rtcm.ssr[i].deph[2]<< 

						mytime[i] = myconv[0]->rtcm.ssr[i].t0[0];
						if (i == 31)  recodfile << endl;
					}
				}
			}
		}

		/*double dclk = myconv[0]->rtcm.ssr[0].dclk[0];
		printf("*****%.4f*****    \n",dclk);*/
		unlock(&mysvr.lock);
		/* get stream server status */
		strsvrstat(&mysvr, strstat, strbyte, bps, strmsg);

		/* show stream server status */
		for (i = 0, p = buff; i < 1; i++)
		{
			p += sprintf(p, "%c", ss[strstat[i] + 1]);
		}

		printf("%s [%s] %10d B %7d bps %s\n",
			time_str(utc2gpst(timeget()), 0), buff, strbyte[0], bps[0], strmsg);

		sleepms(dispint);
	}
}


void testSRIF();

int main(int argc, char* argv[])
{
	g_mylog = fopen("mytrace.log", "w+");

	//double test = fmod(4.0,2.2);

	prcopt_t prcopt = prcopt_default;
	solopt_t solopt = solopt_default;

	prcopt.prn[2] = 0.02 / 3600;  //zwd的大约是2cm每小时

	//set for prcopt
	prcopt.navsys = SYS_GPS;
	prcopt.tidecorr = 2;  // 改正固体潮汐，极移改正,海洋潮汐等
	prcopt.elmin = 5.0*D2R;
	prcopt.posopt[0] = 1; prcopt.posopt[1] = 1; prcopt.posopt[2] = 1; prcopt.posopt[3] = 1; prcopt.posopt[4] = 1; prcopt.posopt[5] = 1;
	prcopt.ionoopt = IONOOPT_IFLC;
	prcopt.tropopt = TROPOPT_EST; //估计ZTD
	prcopt.mode = PMODE_PPP_KINEMA;
	prcopt.modear = 0;
	prcopt.sateph = EPHOPT_PREC;
	prcopt.nf = 1;  //设置单频

	char prog[64] = "                                WHU_PPP                        ";
	//set for solopt
	solopt.datum = 0;
	solopt.posf = SOLF_ENU;
	//solopt.prog = prog;
	solopt.times = TIMES_GPST;
	solopt.sstat = 0;

	int stat = 0;

	//prcopt.exsats[24] = 1;  //排除掉17号卫星
	//prcopt.exsats[26] = 1;  //排除掉13号卫星
	//prcopt.exsats[29] = 1;  //排除掉13号卫星
	//prcopt.exsats[0] = 1;  //排除掉13号卫星
	//prcopt.exsats[2] = 1;  //排除掉13号卫星
	//prcopt.exsats[8] = 1;  //排除掉13号卫星

	//filopt_t filopt={""};
	//gtime_t ts={0},te={0};
	//double ti=0.0,tu=0.0;
	//int i,n=0;
	//char infile_[5][1024]={""},*infile[5],outfile[1024];
	//char *rov,*base,*p,*q,*r;  
	//// set input/output files
	//for (i=0;i<5;i++) 
	//{
	//	infile[i]=infile_[i];
	//}

	//// set input file infile[n]
	//infile[0] = "C:\\Users\\lizhen\\Desktop\\ppp\\ppprt\\demo\\rinex\\wuhn0260.15o";
	//infile[1] = "C:\\Users\\lizhen\\Desktop\\ppp_exp\\brdc0260.15n";
	////infile[2] = "C:\\Users\\lizhen\\Desktop\\ppp_exp\\igs18291.sp3";
	////infile[3] = "C:\\Users\\lizhen\\Desktop\\ppp_exp\\igs18291.clk_30s";
	////infile[4] = "C:\\Users\\lizhen\\Desktop\\ppp_exp\\igsg0260.15i";
	//n = 2;  // 实际输入的文件的个数infile
	//
	//
	//
	//std::string  rov_filePath ="1" ;
	//std::string  base_filePath="0";

	//rov = new char[strlen(rov_filePath.c_str())];  // rov obs file path
	//base = new char[strlen(base_filePath.c_str())]; // base obs file path

	//for ( p=(char*)rov_filePath.c_str(),r=rov;*p;p=q+2) 
	//{
	//	if (!(q=strstr(p,"\r\n"))) 
	//	{
	//		if (*p!='#') strcpy(r,p); break;
	//	}
	//	else if (*p!='#') 
	//	{
	//		strncpy(r,p,q-p); r+=q-p;
	//		strcpy(r++," ");
	//	}
	//}
	//for (p=(char*)base_filePath.c_str(),r=base;*p;p=q+2) 
	//{		
	//	if (!(q=strstr(p,"\r\n")))
	//	{
	//		if (*p!='#') strcpy(r,p); break;
	//	}
	//	else if (*p!='#') 
	//	{
	//		strncpy(r,p,q-p); r+=q-p;
	//		strcpy(r++," ");
	//	}
	//}

	//stat=postpos(ts,te,ti,tu,&prcopt,&solopt,&filopt,infile,n,outfile,rov,base);

	//if(rov != NULL ) { delete[] rov ;  rov = NULL;}
	//if(base != NULL ){ delete[] base;  base = NULL;}

	double myxyz[3] = { -2279833.698, 5004717.650, 3219790.320 };
	double myxyz0[3] = { -2279828.943, 5004706.506, 3219777.437 };
	double pos[3];
	double enu[3] = { 0.0 };
	ecef2pos(myxyz0, pos);
	ecef2enu(pos, myxyz, enu);

	//getEPH(NULL);
	myRTpos(NULL);

	//getBinexData(NULL,0);
	//getUB380Data(NULL,0);
	//getBrdcCorrection( NULL,0);

	//stat =  mypost( &prcopt, &solopt);

	return stat;

}

int getEPH(char* ntrip_str){
	if (ntrip_str == NULL)  { ntrip_str = "001@whu:123|119.97.244.11:2101/GK03"; }  //"gxwang:wang123@58.49.58.149:2101/CLK91"
	int dispint = 1000;
	//for the real time test
	strsvr_t  mysvr;
	int opts[] = { 10000, 10000, 2000, 32768, 10, 0, 30 };
	double stapos[3] = { 0 };

	char *paths[MAXSTR] = { "station:123456|59.175.223.165:2101/NTSCR9" };  //Get observation and navigation data from xian
	paths[0] = ntrip_str;
	paths[1] = "mytest.dat";

	int types[MAXSTR] = { STR_NTRIPCLI, STR_FILE };
	strconv_t *myconv[MAXSTR] = { NULL };
	char *msgout = "1019,1004,1003,1238,1078,1060,1066";
	char* opt = "";
	int strstat[MAXSTR] = { 0 };
	int strbyte[MAXSTR] = { 0 };
	int bps[MAXSTR] = { 0 };
	char strmsg[MAXSTRMSG] = "";
	const char ss[] = { 'E', '-', 'W', 'C', 'C' };
	char buff[256], *p;
	int i = 0, j = 0, n = 1;  // n的值至少是1
	myconv[0] = strconvnew(STRFMT_BINEX, STRFMT_RTCM3, msgout, 0, 0, opt);
	fstream  recodfile("log/igs.rec", ios::out);
	nav_t mynav = { { 0 } };
	obs_t myobs = { { 0 } };
	strsvrinit(&mysvr, n);

	int test_svr = strsvrstart(&mysvr, opts, types, paths, myconv, NULL, stapos);
	if (!test_svr)
	{
		printf("stream server start error\n");
		return -1;
	}

	gtime_t  mytime[32] = { { 0 } };
	char timestr[20] = { 0 }, satstr[3] = { 0 };
	int tag = 0, prn = -1, sys = -1;
	for (intrflg = 0; !intrflg;)
	{
		lock(&mysvr.lock);
		mynav = myconv[0]->raw.nav;
		myobs = myconv[0]->raw.obs;
		printf("星历信息:%d---%d\n", myconv[0]->raw.ephsat, mynav.n);
		for (i = 0; i < mynav.n; i++)
		{
			printf("%d  ", mynav.eph[i].sat);
		}
		for (i = 0; i < mynav.n; i++)
		{
			sys = satsys(mynav.eph[i].sat, &prn);
			if (sys != SYS_GPS)
			{
				continue;
			}
			satno2id(mynav.eph[i].sat, satstr);
			printf("%s:%.4f ", satstr, mynav.eph[i].A);
		}
		printf("\n");

		time2str(myobs.data[0].time, timestr, 3);
		printf("观测值信息:\n");
		printf("epochtime:%s             \n", timestr);
		for (i = 0; i < myobs.n; i++)
		{
			sys = satsys(myobs.data[i].sat, &prn);
			if (sys != SYS_GPS)
			{
				continue;
			}
			satno2id(myobs.data[i].sat, satstr);
			printf("sat: %s obs:%.3f  %.3f %.3f %.3f \n", satstr,
				myobs.data[i].P[0], myobs.data[i].P[1],
				myobs.data[i].L[0], myobs.data[i].L[1]);
		}
		printf("\n");

		unlock(&mysvr.lock);
		/* get stream server status */
		strsvrstat(&mysvr, strstat, strbyte, bps, strmsg);

		/* show stream server status */
		for (i = 0, p = buff; i < 1; i++)
		{
			p += sprintf(p, "%c", ss[strstat[i] + 1]);
		}

		printf("%s [%s] %10d B %7d bps %s\n",
			time_str(utc2gpst(timeget()), 0), buff, strbyte[0], bps[0], strmsg);

		sleepms(dispint);
	}
	return 0;
}

