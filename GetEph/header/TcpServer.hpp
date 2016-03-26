#ifndef _TCPSERVER_H_
#define _TCPSERVER_H_

#include <vector>

#include "rtklib.h"
#include "SafeLog.hpp"
#define MAXSTR 5

class TcpServer
{
public:
	TcpServer()
	{
		io_types[0] = STR_NONE;
		io_types[1] = STR_TCPSVR;
		n_eph = NSATGPS + NSATGAL + NSATCMP + NSATGLO;
		n_geph = NSATGLO;
		opts = new int[7]{10000, 10000, 2000, 32768, 10, 0, 30};
		nout = 1;
		strsvrinit(&mysvr, nout);
		initConv();
	}

	~TcpServer()
	{
	}

	void setIOpaths(char* i_path, char* o_path)
	{
		io_paths[0] = "";
		io_paths[1] = o_path;
	}

	int start()
	{
		int result = 0;
		result = strsvrstart(&mysvr, opts, io_types, io_paths, &out_conv, NULL, NULL);
		if (result == 1)
		{
			cout << port << "打开成功" << endl;
		}
		else
		{
			cout << port << "打开失败" << endl;
		}
		return result;
	}
	void close()
	{
		strsvrstop(&mysvr, NULL);
		this->~TcpServer();
	}

	void setPort(int port)
	{
		this->port = port;
		memset(cport, 0, 6 * sizeof(char));
		sprintf(cport, ":%d%c", port, '\0');
		setIOpaths(NULL, cport);
	}

	void push_data(rtcm_t data)
	{
		lock(&mysvr.lock);

		memcpy(out_conv->out.nav.eph, data.nav.eph, sizeof(data.nav.eph)*n_eph);
		memcpy(out_conv->out.nav.geph, data.nav.geph, sizeof(data.nav.geph)*n_geph);

		unlock(&mysvr.lock);
	}

private:
	strsvr_t  mysvr;

	char *io_paths[2];
	int io_types[2];
	int *opts;
	char *msgout;

	int itype, otype;
	char* opt;

	int port;
	char cport[6];

	strconv_t* out_conv;

	void initConv()
	{
		out_conv = strconvnew(STRFMT_RTCM3, STRFMT_RTCM3, "1019(2),1020(2),1045(2),1047(2)", 0, 0, "");

		out_conv->out.nav.eph = (eph_t*)malloc(sizeof(eph_t)* n_eph);
		out_conv->out.nav.geph = (geph_t *)malloc(sizeof(geph_t)* n_geph);

		memset(out_conv->out.nav.eph, 1, sizeof(eph_t)*n_eph);
		for (int i = 0; i < n_eph; i++)
		{
			out_conv->out.nav.eph[i].sat = i + 1;
		}
		memset(out_conv->out.nav.geph, 1, sizeof(geph_t)*n_geph);
		for (int i = 0; i < n_geph; i++)
		{
			out_conv->out.nav.geph[i].sat = NSATGPS + i + 1;
		}
	}

	int nout;

	int n_eph, n_geph;

};

#endif
