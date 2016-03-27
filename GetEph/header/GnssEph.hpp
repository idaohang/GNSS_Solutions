#ifndef _GNSSEPH_H_
#define _GNSSEPH_H_

#include <vector>

#include "rtklib.h"
#include "SafeLog.hpp"

#define MAXSTR 5
class GNSSEPH
{
public:
	GNSSEPH()
	{
		n_eph = NSATGPS + NSATGAL + NSATCMP + NSATGLO;
		n_geph = NSATGLO;

		gnss_sat_num = n_eph;

		initRtcm();
	}

	~GNSSEPH()
	{

	}

	void addEph(int sat, const eph_t* eph) //gps/gal/cmp eph
	{
		memcpy(&rtcm_data.nav.eph[sat - 1], eph, sizeof(eph_t));
	}

	void addGeph(int sat, const geph_t* geph)
	{
		int gsat = sat - NSATGPS;
		memcpy(&rtcm_data.nav.geph[gsat - 1], geph, sizeof(geph_t));
	}

	void log()
	{
		int sat = 0;
		SafeLog::logMsg("All the eph:");

		for (int i = 0; i < gnss_sat_num; i++)
		{

			if (32 <= i&&i <= 55)
			{
				sat = rtcm_data.nav.geph[i - NSATGPS].sat;
			}
			else
			{
				sat = rtcm_data.nav.eph[i].sat;
			}

			SafeLog::logInt(sat);
		}

		SafeLog::logEndl();
	}

	void print()
	{
		int sat = 0;
		cout << "All the eph:" << endl;

		for (int i = 0; i < gnss_sat_num; i++)
		{

			if (32 <= i&&i <= 55)
			{
				sat = rtcm_data.nav.geph[i - NSATGPS].sat;
			}
			else
			{
				sat = rtcm_data.nav.eph[i].sat;
			}

			if (sat != 0)
			{
				cout << sat << "  ";
			}
		}

		cout << endl;
	}

	strconv_t* myconv[MAXSTR];

	void setTime(gtime_t time)
	{
		this->time = time;
		rtcm_data.time = time;
	}

	rtcm_t getRtcm_data()
	{
		return rtcm_data;
	}
	
private:
	int gnss_sat_num;
	gtime_t time;
	int n_eph, n_geph;
	rtcm_t rtcm_data;
	void initRtcm()
	{
		rtcm_data.nav.eph = (eph_t*)malloc(sizeof(eph_t)* n_eph);
		rtcm_data.nav.geph = (geph_t *)malloc(sizeof(geph_t)* n_geph);

		memset(rtcm_data.nav.eph, 0, sizeof(eph_t)*n_eph);
		memset(rtcm_data.nav.geph, 0, sizeof(geph_t)*n_geph);
	}
};


#endif
