#include <iostream>
#include <string>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <vector>
#include "../header/StrSvr.hpp"
#include "../header/GnssEph.hpp"
#include "../header/rtklib.h"
#include "../header/TcpServer.hpp"
#include "../header/SafeLog.hpp"

using namespace std;

FILE* g_mylog = NULL;

#define MAXSTR      5 

GNSSEPH gnssEph;

TcpServer server;

thread_t getEphThread;

#ifdef WIN32
static DWORD WINAPI  getEph(void *arg)
#else
static void *getEph(void *arg)
#endif
{
	rtcm_t data;
	int dispint = 1000; //设置间隔1000ms
	int t1;
	while (true)
	{
		t1 = tickget();

		for (int i=0;i<StrSvr::getStrSvrs().size();i++)
		{
			StrSvr strsvr=StrSvr::getStrSvrs().at(i);

			data = *strsvr.getRtcm();

			gnssEph.setTime(data.time);

			for (int i = 0; i < MAXSAT; i++)
			{
				int sat = 0;
				if (32 <= i && i <= 55)
				{
					if ((sat = data.nav.geph[i - NSATGPS].sat) != 0) // glonass eph
					{
						gnssEph.addGeph(sat, &data.nav.geph[sat - NSATGPS - 1]);
						continue;
					}
				}
				else
				{
					if ((sat = data.nav.eph[i].sat) != 0)
					{
						gnssEph.addEph(sat, &data.nav.eph[sat - 1]); // gps/galileo/cmp eph
						continue;
					}
				}
			}
		}
		gnssEph.log();

		server.push_data(gnssEph.getRtcm_data());

		sleepms(dispint - (tickget() - t1));
	}

	return 0;
}

void printHelp()
{
	cout << "q:退出" << endl;
	cout << "h:帮助" << endl;
	cout << "p:打印" << endl;
}

void timeDiff(time_t start, time_t end, char* out)
{
	double dt = end - start;
	double ds = 3600 * 24;
	double hs = 3600;
	double ms = 60;
	int day = (int)(dt / ds);
	int hour = (int)((dt - day*ds) / hs);
	int minute = (int)((dt - day*ds - hour*hs) / ms);
	double second = dt - day*ds - hour*hs - minute*ms;

	sprintf(out, "%dd %dh %dm %.1lfs", day, hour, minute, second);
}

void startGetEph()
{
#ifdef WIN32
	if (!(getEphThread = CreateThread(NULL, 0, getEph, NULL, 0, NULL))) {
#else
	if (pthread_create(&getEphThread, NULL, getEph, NULL)) {
#endif
		cout << "创建线程失败" << endl;
	}
}

void stopGetEph()
{
#ifdef WIN32
	WaitForSingleObject(getEphThread, 1000);
	CloseHandle(getEphThread);
#else
	pthread_join(getEphThread, NULL);
#endif
}

int main(int argc, char* argv[])
{
#ifdef WIN32
	system("title get eph from caster mountpoints");
#else

#endif
	system("mkdir log");
	g_mylog = fopen("log/mytrace.log", "w+");

	time_t startTime, endTime;
	time(&startTime);

	int port = 0;
	cout << "输入端口号：";//输入端口号
	cin >> port;

	server.setPort(port);

	StrSvr::loadFromFile("sourceTable.txt");

	if (!server.start())
	{
		StrSvr::closeStrSvrs();
		return 0;
	}

	StrSvr::startStrSves();
	startGetEph();

	char cmd = NULL;
	printHelp();
	while (true)
	{
		cmd = getchar();
		switch (cmd)
		{
		case 'q':
		case 'Q':
			cout << "waitting...";
			stopGetEph();
			StrSvr::closeStrSvrs();
			server.close();
			return 0;
			break;
		case 'h':
		case 'H':
			printHelp();
			break;
		case 'p':
		case 'P':
			time(&endTime);
			char runtime[20] = { 0 };
			timeDiff(startTime, endTime, runtime);
			cout << "运行时间：" << runtime << endl;
			cout << "Tcp Server 输出端口：" << port << endl;
			gnssEph.print();
			break;
		}
	}
	return 0;
}
