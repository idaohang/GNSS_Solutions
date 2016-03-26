#ifndef _SAFELOG_H_
#define _SAFELOG_H_

#include <iostream>
#include <time.h>
#include <stdio.h>

using namespace std;

class SafeLog
{
public :

	SafeLog()
	{
	}

	~SafeLog()
	{
	}

	static void logMsg(char* msg)
	{
		logTime();
		fprintf(getfp(), "%s\n", msg);
	}
	static void logInt(int n)
	{
		fprintf(getfp(), "%d  ", n);
	}
	static void logEndl()
	{
		fprintf(getfp(), "\n");
	}
	static void logTime()
	{
		char t[64] = { 0 };
		getTime(t);
		fprintf(getfp(), "%s : ", t);
	}
	static void close()
	{
		fclose(getfp());
	}

private:
	
	static void getTime(char* ct)
	{
		time_t t = time(0);
		char ctime[64];
		strftime(ctime, sizeof(ctime), "%Y/%m/%d %H:%M:%S", localtime(&t));
		memcpy(ct, ctime, 20 * sizeof(char));
	}

	static FILE* getfp()
	{
		static FILE* logfp = fopen("log/safelog.log", "a");
		fflush(logfp);
		return logfp;
	}
};

#endif
