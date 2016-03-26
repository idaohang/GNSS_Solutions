#ifndef _STRSVR_H_
#define _STRSVR_H_

#include <vector>

#include "rtklib.h"
#include "SafeLog.hpp"

using namespace std;
#define MAXSTR      5                  /* max number of streams */

class StrSvr
{
public:
	StrSvr()
	{
		defaultSet();
	}
	~StrSvr()
	{
	}

	int start()
	{
		int result = 0;
		char msg[20] = { 0 };
		myconv[0] = strconvnew(itype, otype, msgout, 0, 0, opt);
		strsvrinit(&mysvr, nout);

		result = strsvrstart(&mysvr, opts, types, paths, myconv, NULL, NULL);
		if (result == 0)
		{
			cout << getName() << "打开失败" << endl;
			sprintf(msg, "%s:打开失败。", getName());
		}
		else
		{
			sprintf(msg, "%s:打开成功。", getName());
		}

		SafeLog::logMsg(msg);
		return result;
	}
	
	void close()
	{
		strsvrstop(&mysvr, NULL);
		this->~StrSvr();
	}

	const rtcm_t* getRtcm()
	{
		rtcm_t* out;
		lock(&mysvr.lock);
		out = &(myconv[0]->out);
		unlock(&mysvr.lock);
		log();
		return out;
	}
	
	void setInputPath(char* path)
	{
		paths[0] = (char*)malloc(sizeof(char)* 100);
		memcpy(paths[0], path, 100 * sizeof(char));

		int count = 0;
		char *p = paths[0];

		for (count = 0; p + count != ""; count++)
		{
			if (*(p + count) == '/')
			{
				name = paths[0] + count + 1;
				break;
			}
		}
	}
	
	void setOutputPaths(char* path1, char* path2, char* path3, int n)
	{
		paths[1] = path1;
		paths[2] = path2;
		paths[3] = path3;
		nout = n;
	}

	const char* getName()
	{
		return name;
	}
	
	void setIOtype(int itype, int otype)
	{
		this->itype = itype;
		this->otype = otype;
	}
	
	void setStrIOtype(int Itype, int Otype1, int Otype2, int Otype3)
	{
		memset(types, 0, sizeof(int)*MAXSTR);
		types[0] = Itype;
		types[1] = Otype1;
		types[2] = Otype2;
		types[3] = Otype3;
	}

	static int loadFromFile(char* fileDir)
	{
		ifstream fin(fileDir);
		char str_dir[100] = { 0 }, str_type[10] = { 0 };
		int type = 0;

		while (!fin.eof())
		{
			fin >> str_dir >> str_type;

			if (strncmp(str_type, "RTCM3", 5) == 0)
			{
				type = STRFMT_RTCM3;
			}
			else if (strncmp(str_type, "BINEX", 5) == 0)
			{
				type = STRFMT_BINEX;
			}

			StrSvr str = StrSvr::getStrSvr(str_dir, type);

			getStrSvrs().push_back(str);
		}

		return getStrSvrs().size();
	}
	
	static StrSvr getStrSvr(char* dir, int type)
	{
		StrSvr strsvr;
		strsvr.setInputPath(dir);
		strsvr.setIOtype(type, STRFMT_RTCM3);
		return strsvr;
	}
	
	static vector<StrSvr>& getStrSvrs()
	{
		static vector<StrSvr> strsvrs;
		return strsvrs;
	}

	static void addStrSvrs(StrSvr strsvr)
	{
		getStrSvrs().push_back(strsvr);
	}
	
	static void startStrSves()
	{
		for (int n = 0; n < getStrSvrs().size(); ++n)
		{
			getStrSvrs().at(n).start();
		}
	}
	
	static void closeStrSvrs()
	{
		for (int n = 0; n < getStrSvrs().size(); ++n)
		{
			getStrSvrs().at(n).close();
		}
	}

private:
	char* name;
	char *paths[MAXSTR];
	int types[MAXSTR]; //输入输出流类型
	int *opts;  // stream options （timeout。。。。）
	char* opt;
	char *msgout;
	int nout; // 输出流个数
	int itype, otype; // 输入输出数据类型

	strsvr_t  mysvr;
	strconv_t *myconv[MAXSTR];

	void defaultSet()
	{
		types[0] = STR_NTRIPCLI;
		types[1] = STR_NONE;

		msgout = "1004(1),1012(1),1006(10),1008(10),1019(15),1033(10),1045(15),1046(15),1077(1),1087(1),1097(1),1107(1),1127(1)";
		opts = new int[7]{ 10000, 10000, 2000, 32768, 10, 0, 30 };
		itype = STRFMT_RTCM3;
		otype = STRFMT_RTCM3;
		opt = "";
		nout = 1;
		setOutputPaths("mytest", NULL, NULL, nout);
	}

	void log()
	{
		int sat = 0;

		char msg[20] = { 0 };
		sprintf(msg, "%s::%d", getName(), myconv[0]->out.obs.n);
		SafeLog::logMsg(msg);
		//cout << msg << endl;

		for (int i = 0; i < myconv[0]->out.nav.n; i++)
		{
			if (32 <= i&&i <= 55)
			{
				if ((sat = myconv[0]->out.nav.geph[i - NSATGPS].sat) != 0)
				{
					SafeLog::logInt(sat);
				}
			}
			else
			{
				if ((sat = myconv[0]->out.nav.eph[i].sat) != 0)
				{
					SafeLog::logInt(sat);
				}
			}
			if (sat != 0)
			{
				//cout << sat << "  ";
			}
		}

		//cout << endl;
		SafeLog::logEndl();
	}
};
#endif
