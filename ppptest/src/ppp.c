/*------------------------------------------------------------------------------
* ppp.c : precise point positioning
*
*          Copyright (C) 2010-2013 by T.TAKASU, All rights reserved.
*
* options : -DIERS_MODEL use IERS tide model
*
* references :
*     [1] D.D.McCarthy, IERS Technical Note 21, IERS Conventions 1996, July 1996
*     [2] D.D.McCarthy and G.Petit, IERS Technical Note 32, IERS Conventions
*         2003, November 2003
*     [3] D.A.Vallado, Fundamentals of Astrodynamics and Applications 2nd ed,
*         Space Technology Library, 2004
*     [4] J.Kouba, A Guide to using International GNSS Service (IGS) products,
*         May 2009
*     [5] RTCM Paper, April 12, 2010, Proposed SSR Messages for SV Orbit Clock,
*         Code Biases, URA
*     [6] MacMillan et al., Atmospheric gradients and the VLBI terrestrial and
*         celestial reference frames, Geophys. Res. Let., 1997
*     [7] G.Petit and B.Luzum (eds), IERS Technical Note No. 36, IERS
*         Conventions (2010), 2010
*
* version : $Revision:$ $Date:$
* history : 2010/07/20 1.0  new
*                           added api:
*                               tidedisp()
*           2010/12/11 1.1  enable exclusion of eclipsing satellite
*           2012/02/01 1.2  add gps-glonass h/w bias correction
*                           move windupcorr() to rtkcmn.c
*           2013/03/11 1.3  add otl and pole tides corrections
*                           involve iers model with -DIERS_MODEL
*                           change initial variances
*                           suppress acos domain error
*           2013/09/01 1.4  pole tide model by iers 2010
*                           add mode of ionosphere model off
*           2014/05/23 1.5  add output of trop gradient in solution status
*           2014/10/13 1.6  fix bug on P0(a[3]) computation in tide_oload()
*                           fix bug on m2 computation in tide_pole()
*-----------------------------------------------------------------------------*/
#include "rtklib.h"

static const char rcsid[] = "$Id:$";

#define SQR(x)      ((x)*(x))
#define MIN(x,y)    ((x)<=(y)?(x):(y))

#define AS2R        (D2R/3600.0)    /* arc sec to radian */
#define GME         3.986004415E+14 /* earth gravitational constant */
#define GMS         1.327124E+20    /* sun gravitational constant */
#define GMM         4.902801E+12    /* moon gravitational constant */

/* initial variances */
#define VAR_POS     SQR(100.0)      /*   receiver position (m^2) */
#define VAR_CLK     SQR(100.0)      /*   receiver clock (m^2) */
#define VAR_ZTD     SQR(  0.3)      /*   ztd (m^2) */
#define VAR_GRA     SQR(0.001)      /*   gradient (m^2) */
#define VAR_BIAS    SQR(100.0)      /*   phase-bias (m^2) */

#define VAR_IONO_OFF SQR(10.0)      /* variance of iono-model-off */

#define ERR_SAAS    0.3             /* saastamoinen model error std (m) */
#define ERR_BRDCI   0.5             /* broadcast iono model error factor */
#define ERR_CBIAS   0.3             /* code bias error std (m) */
#define REL_HUMI    0.7             /* relative humidity for saastamoinen model */

#define NP(opt)     ((opt)->dynamics?9:3) /* number of pos solution */

#define IC(s,opt)   (NP(opt)+(s))      /* state index of clocks (s=0:gps,1:glo:,2:bds) */
#define IT(opt)     (IC(0,opt)+NSYS)   /* state index of tropos */

#define NR(opt)     (IT(opt)+((opt)->tropopt<TROPOPT_EST?0:((opt)->tropopt==TROPOPT_EST?1:3)))
/* number of solutions */
#define IB(s,opt)   (NR(opt)+(s)-1)    /* state index of phase bias */
#define NX(opt)     (IB(MAXSAT,opt)+1) /* number of estimated states */



/* function prototypes -------------------------------------------------------*/
#ifdef IERS_MODEL
extern int dehanttideinel_(double *xsta, int *year, int *mon, int *day,
	double *fhr, double *xsun, double *xmon,
	double *dxtide);
#endif

/* output solution status for PPP --------------------------------------------*/
extern void pppoutsolstat(rtk_t *rtk, int level, FILE *fp)
{
	ssat_t *ssat;
	double tow, pos[3], vel[3], acc[3];
	int i, j, week, nfreq = 1;
	char id[32];

	if (level <= 0 || !fp) return;

	trace(3, "pppoutsolstat:\n");

	tow = time2gpst(rtk->sol.time, &week);

	/* receiver position */
	fprintf(fp, "$POS,%d,%.3f,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n", week, tow,
		rtk->sol.stat, rtk->x[0], rtk->x[1], rtk->x[2], 0.0, 0.0, 0.0);

	/* receiver velocity and acceleration */
	if (rtk->opt.dynamics) {
		ecef2pos(rtk->sol.rr, pos);
		ecef2enu(pos, rtk->x + 3, vel);
		ecef2enu(pos, rtk->x + 6, acc);
		fprintf(fp, "$VELACC,%d,%.3f,%d,%.4f,%.4f,%.4f,%.5f,%.5f,%.5f,%.4f,%.4f,%.4f,%.5f,%.5f,%.5f\n",
			week, tow, rtk->sol.stat, vel[0], vel[1], vel[2], acc[0], acc[1], acc[2],
			0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
	}
	/* receiver clocks */
	i = IC(0, &rtk->opt);
	fprintf(fp, "$CLK,%d,%.3f,%d,%d,%.3f,%.3f,%.3f,%.3f\n",
		week, tow, rtk->sol.stat, 1, rtk->x[i] * 1E9 / CLIGHT, rtk->x[i + 1] * 1E9 / CLIGHT,
		0.0, 0.0);

	/* tropospheric parameters */
	if (rtk->opt.tropopt == TROPOPT_EST || rtk->opt.tropopt == TROPOPT_ESTG) {
		i = IT(&rtk->opt);
		fprintf(fp, "$TROP,%d,%.3f,%d,%d,%.4f,%.4f\n", week, tow, rtk->sol.stat,
			1, rtk->x[i], 0.0);
	}
	if (rtk->opt.tropopt == TROPOPT_ESTG) {
		i = IT(&rtk->opt);
		fprintf(fp, "$TRPG,%d,%.3f,%d,%d,%.5f,%.5f,%.5f,%.5f\n", week, tow,
			rtk->sol.stat, 1, rtk->x[i + 1], rtk->x[i + 2], 0.0, 0.0);
	}
	if (rtk->sol.stat == SOLQ_NONE || level <= 1) return;

	/* residuals and status */
	for (i = 0; i < MAXSAT; i++) {
		ssat = rtk->ssat + i;
		if (!ssat->vs) continue;
		satno2id(i + 1, id);
		for (j = 0; j < nfreq; j++) {
			fprintf(fp, "$SAT,%d,%.3f,%s,%d,%.1f,%.1f,%.4f,%.4f,%d,%.0f,%d,%d,%d,%d,%d,%d\n",
				week, tow, id, j + 1, ssat->azel[0] * R2D, ssat->azel[1] * R2D,
				ssat->resp[j], ssat->resc[j], ssat->vsat[j], ssat->snr[j] * 0.25,
				ssat->fix[j], ssat->slip[j] & 3, ssat->lock[j], ssat->outc[j],
				ssat->slipc[j], ssat->rejc[j]);
		}
	}
}
/* solar/lunar tides (ref [2] 7) ---------------------------------------------*/
static void tide_pl(const double *eu, const double *rp, double GMp,
	const double *pos, double *dr)
{
	const double H3 = 0.292, L3 = 0.015;
	double r, ep[3], latp, lonp, p, K2, K3, a, H2, L2, dp, du, cosp, sinl, cosl;
	int i;

	trace(4, "tide_pl : pos=%.3f %.3f\n", pos[0] * R2D, pos[1] * R2D);

	if ((r = norm(rp, 3)) <= 0.0) return;

	for (i = 0; i < 3; i++) ep[i] = rp[i] / r;

	K2 = GMp / GME*SQR(RE_WGS84)*SQR(RE_WGS84) / (r*r*r);
	K3 = K2*RE_WGS84 / r;
	latp = asin(ep[2]); lonp = atan2(ep[1], ep[0]);
	cosp = cos(latp); sinl = sin(pos[0]); cosl = cos(pos[0]);

	/* step1 in phase (degree 2) */
	p = (3.0*sinl*sinl - 1.0) / 2.0;
	H2 = 0.6078 - 0.0006*p;
	L2 = 0.0847 + 0.0002*p;
	a = dot(ep, eu, 3);
	dp = K2*3.0*L2*a;
	du = K2*(H2*(1.5*a*a - 0.5) - 3.0*L2*a*a);

	/* step1 in phase (degree 3) */
	dp += K3*L3*(7.5*a*a - 1.5);
	du += K3*(H3*(2.5*a*a*a - 1.5*a) - L3*(7.5*a*a - 1.5)*a);

	/* step1 out-of-phase (only radial) */
	du += 3.0 / 4.0*0.0025*K2*sin(2.0*latp)*sin(2.0*pos[0])*sin(pos[1] - lonp);
	du += 3.0 / 4.0*0.0022*K2*cosp*cosp*cosl*cosl*sin(2.0*(pos[1] - lonp));

	dr[0] = dp*ep[0] + du*eu[0];
	dr[1] = dp*ep[1] + du*eu[1];
	dr[2] = dp*ep[2] + du*eu[2];

	trace(5, "tide_pl : dr=%.3f %.3f %.3f\n", dr[0], dr[1], dr[2]);
}
/* displacement by solid earth tide (ref [2] 7) ------------------------------*/
static void tide_solid(const double *rsun, const double *rmoon,
	const double *pos, const double *E, double gmst, int opt,
	double *dr)
{
	double dr1[3], dr2[3], eu[3], du, dn, sinl, sin2l;

	trace(3, "tide_solid: pos=%.3f %.3f opt=%d\n", pos[0] * R2D, pos[1] * R2D, opt);

	/* step1: time domain */
	eu[0] = E[2]; eu[1] = E[5]; eu[2] = E[8];
	tide_pl(eu, rsun, GMS, pos, dr1);
	tide_pl(eu, rmoon, GMM, pos, dr2);

	/* step2: frequency domain, only K1 radial */
	sin2l = sin(2.0*pos[0]);
	du = -0.012*sin2l*sin(gmst + pos[1]);

	dr[0] = dr1[0] + dr2[0] + du*E[2];
	dr[1] = dr1[1] + dr2[1] + du*E[5];
	dr[2] = dr1[2] + dr2[2] + du*E[8];

	/* eliminate permanent deformation */
	if (opt & 8) {
		sinl = sin(pos[0]);
		du = 0.1196*(1.5*sinl*sinl - 0.5);
		dn = 0.0247*sin2l;
		dr[0] += du*E[2] + dn*E[1];
		dr[1] += du*E[5] + dn*E[4];
		dr[2] += du*E[8] + dn*E[7];
	}
	trace(5, "tide_solid: dr=%.3f %.3f %.3f\n", dr[0], dr[1], dr[2]);
}
/* displacement by ocean tide loading (ref [2] 7) ----------------------------*/
static void tide_oload(gtime_t tut, const double *odisp, double *denu)
{
	const double args[][5] = {
		{ 1.40519E-4, 2.0, -2.0, 0.0, 0.00 },  /* M2 */
		{ 1.45444E-4, 0.0, 0.0, 0.0, 0.00 },  /* S2 */
		{ 1.37880E-4, 2.0, -3.0, 1.0, 0.00 },  /* N2 */
		{ 1.45842E-4, 2.0, 0.0, 0.0, 0.00 },  /* K2 */
		{ 0.72921E-4, 1.0, 0.0, 0.0, 0.25 },  /* K1 */
		{ 0.67598E-4, 1.0, -2.0, 0.0, -0.25 },  /* O1 */
		{ 0.72523E-4, -1.0, 0.0, 0.0, -0.25 },  /* P1 */
		{ 0.64959E-4, 1.0, -3.0, 1.0, -0.25 },  /* Q1 */
		{ 0.53234E-5, 0.0, 2.0, 0.0, 0.00 },  /* Mf */
		{ 0.26392E-5, 0.0, 1.0, -1.0, 0.00 },  /* Mm */
		{ 0.03982E-5, 2.0, 0.0, 0.0, 0.00 }   /* Ssa */
	};
	const double ep1975[] = { 1975, 1, 1, 0, 0, 0 };
	double ep[6], fday, days, t, t2, t3, a[5], ang, dp[3] = { 0 };
	int i, j;

	trace(3, "tide_oload:\n");

	/* angular argument: see subroutine arg.f for reference [1] */
	time2epoch(tut, ep);
	fday = ep[3] * 3600.0 + ep[4] * 60.0 + ep[5];
	ep[3] = ep[4] = ep[5] = 0.0;
	days = timediff(epoch2time(ep), epoch2time(ep1975)) / 86400.0;
	t = (27392.500528 + 1.000000035*days) / 36525.0;
	t2 = t*t; t3 = t2*t;

	a[0] = fday;
	a[1] = (279.69668 + 36000.768930485*t + 3.03E-4*t2)*D2R; /* H0 */
	a[2] = (270.434358 + 481267.88314137*t - 0.001133*t2 + 1.9E-6*t3)*D2R; /* S0 */
	a[3] = (334.329653 + 4069.0340329577*t - 0.010325*t2 - 1.2E-5*t3)*D2R; /* P0 */
	a[4] = 2.0*PI;

	/* displacements by 11 constituents */
	for (i = 0; i < 11; i++) {
		ang = 0.0;
		for (j = 0; j < 5; j++) ang += a[j] * args[i][j];
		for (j = 0; j < 3; j++) dp[j] += odisp[j + i * 6] * cos(ang - odisp[j + 3 + i * 6] * D2R);
	}
	denu[0] = -dp[1];
	denu[1] = -dp[2];
	denu[2] = dp[0];

	trace(5, "tide_oload: denu=%.3f %.3f %.3f\n", denu[0], denu[1], denu[2]);
}
/* iers mean pole (ref [7] eq.7.25) ------------------------------------------*/
static void iers_mean_pole(gtime_t tut, double *xp_bar, double *yp_bar)
{
	const double ep2000[] = { 2000, 1, 1, 0, 0, 0 };
	double y, y2, y3;

	y = timediff(tut, epoch2time(ep2000)) / 86400.0 / 365.25;

	if (y < 3653.0 / 365.25) { /* until 2010.0 */
		y2 = y*y; y3 = y2*y;
		*xp_bar = 55.974 + 1.8243*y + 0.18413*y2 + 0.007024*y3; /* (mas) */
		*yp_bar = 346.346 + 1.7896*y - 0.10729*y2 - 0.000908*y3;
	}
	else { /* after 2010.0 */
		*xp_bar = 23.513 + 7.6141*y; /* (mas) */
		*yp_bar = 358.891 - 0.6287*y;
	}
}
/* displacement by pole tide (ref [7] eq.7.26) --------------------------------*/
static void tide_pole(gtime_t tut, const double *pos, const double *erpv,
	double *denu)
{
	double xp_bar, yp_bar, m1, m2, cosl, sinl;

	trace(3, "tide_pole: pos=%.3f %.3f\n", pos[0] * R2D, pos[1] * R2D);

	/* iers mean pole (mas) */
	iers_mean_pole(tut, &xp_bar, &yp_bar);

	/* ref [7] eq.7.24 */
	m1 = erpv[0] / AS2R - xp_bar*1E-3; /* (as) */
	m2 = -erpv[1] / AS2R + yp_bar*1E-3;

	/* sin(2*theta) = sin(2*phi), cos(2*theta)=-cos(2*phi) */
	cosl = cos(pos[1]);
	sinl = sin(pos[1]);
	denu[0] = 9E-3*sin(pos[0])    *(m1*sinl - m2*cosl); /* de= Slambda (m) */
	denu[1] = -9E-3*cos(2.0*pos[0])*(m1*cosl + m2*sinl); /* dn=-Stheta  (m) */
	denu[2] = -33E-3*sin(2.0*pos[0])*(m1*cosl + m2*sinl); /* du= Sr      (m) */

	trace(5, "tide_pole : denu=%.3f %.3f %.3f\n", denu[0], denu[1], denu[2]);
}
/* tidal displacement ----------------------------------------------------------
* displacements by earth tides
* args   : gtime_t tutc     I   time in utc
*          double *rr       I   site position (ecef) (m)
*          int    opt       I   options (or of the followings)
*                                 1: solid earth tide
*                                 2: ocean tide loading
*                                 4: pole tide
*                                 8: elimate permanent deformation
*          double *erp      I   earth rotation parameters (NULL: not used)
*          double *odisp    I   ocean loading parameters  (NULL: not used)
*                                 odisp[0+i*6]: consituent i amplitude radial(m)
*                                 odisp[1+i*6]: consituent i amplitude west  (m)
*                                 odisp[2+i*6]: consituent i amplitude south (m)
*                                 odisp[3+i*6]: consituent i phase radial  (deg)
*                                 odisp[4+i*6]: consituent i phase west    (deg)
*                                 odisp[5+i*6]: consituent i phase south   (deg)
*                                (i=0:M2,1:S2,2:N2,3:K2,4:K1,5:O1,6:P1,7:Q1,
*                                   8:Mf,9:Mm,10:Ssa)
*          double *dr       O   displacement by earth tides (ecef) (m)
* return : none
* notes  : see ref [1], [2] chap 7
*          see ref [4] 5.2.1, 5.2.2, 5.2.3
*          ver.2.4.0 does not use ocean loading and pole tide corrections
*-----------------------------------------------------------------------------*/
extern void tidedisp(gtime_t tutc, const double *rr, int opt, const erp_t *erp,
	const double *odisp, double *dr)
{
	gtime_t tut;
	double pos[2], E[9], drt[3], denu[3], rs[3], rm[3], gmst, erpv[5] = { 0 };
	int i;
#ifdef IERS_MODEL
	double ep[6],fhr;
	int year,mon,day;
#endif

	trace(3, "tidedisp: tutc=%s\n", time_str(tutc, 0));

	if (erp) geterp(erp, tutc, erpv);

	tut = timeadd(tutc, erpv[2]);

	dr[0] = dr[1] = dr[2] = 0.0;

	if (norm(rr, 3) <= 0.0) return;

	pos[0] = asin(rr[2] / norm(rr, 3));
	pos[1] = atan2(rr[1], rr[0]);
	xyz2enu(pos, E);

	if (opt & 1) { /* solid earth tides */

		/* sun and moon position in ecef */
		sunmoonpos(tutc, erpv, rs, rm, &gmst);

#ifdef IERS_MODEL
		time2epoch(tutc,ep);
		year=(int)ep[0];
		mon =(int)ep[1];
		day =(int)ep[2];
		fhr =ep[3]+ep[4]/60.0+ep[5]/3600.0;

		/* call DEHANTTIDEINEL */
		dehanttideinel_((double *)rr,&year,&mon,&day,&fhr,rs,rm,drt);
#else
		tide_solid(rs, rm, pos, E, gmst, opt, drt);
#endif
		for (i = 0; i < 3; i++) dr[i] += drt[i];
	}
	if ((opt & 2) && odisp) { /* ocean tide loading */
		tide_oload(tut, odisp, denu);
		matmul("TN", 3, 1, 3, 1.0, E, denu, 0.0, drt);
		for (i = 0; i < 3; i++) dr[i] += drt[i];
	}
	if ((opt & 4) && erp) { /* pole tide */
		tide_pole(tut, pos, erpv, denu);
		matmul("TN", 3, 1, 3, 1.0, E, denu, 0.0, drt);
		for (i = 0; i < 3; i++) dr[i] += drt[i];
	}
	trace(5, "tidedisp: dr=%.3f %.3f %.3f\n", dr[0], dr[1], dr[2]);
}
/* exclude meas of eclipsing satellite (block IIA) ---------------------------*/
static void testeclipse(const obsd_t *obs, int n, const nav_t *nav, double *rs)
{
	double rsun[3], esun[3], r, ang, erpv[5] = { 0 }, cosa;
	int i, j;
	const char *type;

	trace(3, "testeclipse:\n");

	/* unit vector of sun direction (ecef) */
	sunmoonpos(gpst2utc(obs[0].time), erpv, rsun, NULL, NULL);
	normv3(rsun, esun);

	for (i = 0; i < n; i++) {
		type = nav->pcvs[obs[i].sat - 1].type;

		if ((r = norm(rs + i * 6, 3)) <= 0.0) continue;
#if 1
		/* only block IIA */
		if (*type&&!strstr(type, "BLOCK IIA")) continue;
#endif
		/* sun-earth-satellite angle */
		cosa = dot(rs + i * 6, esun, 3) / r;
		cosa = cosa<-1.0 ? -1.0 : (cosa>1.0 ? 1.0 : cosa);
		ang = acos(cosa);

		/* test eclipse */
		if (ang<PI / 2.0 || r*sin(ang)>RE_WGS84) continue;

		trace(2, "eclipsing sat excluded %s sat=%2d\n", time_str(obs[0].time, 0),
			obs[i].sat);

		for (j = 0; j < 3; j++) rs[j + i * 6] = 0.0;
	}
}
/* measurement error variance ------------------------------------------------*/
static double varerr(int sat, int sys, double el, int type, const prcopt_t *opt)
{
	double a, b, a2, b2, fact = 1.0;
	double sinel = sin(el);
	int i = sys == SYS_GLO ? 1 : (sys == SYS_GAL ? 2 : 0);

	/* extended error model */
	if (type == 1 && opt->exterr.ena[0])
	{ /* code */
		a = opt->exterr.cerr[i][0];
		b = opt->exterr.cerr[i][1];
		if (opt->ionoopt == IONOOPT_IFLC)
		{
			a2 = opt->exterr.cerr[i][2];
			b2 = opt->exterr.cerr[i][3];
			a = sqrt(SQR(2.55)*a*a + SQR(1.55)*a2*a2);
			b = sqrt(SQR(2.55)*b*b + SQR(1.55)*b2*b2);
		}
	}
	else if (type == 0 && opt->exterr.ena[1])
	{ /* phase */
		a = opt->exterr.perr[i][0];
		b = opt->exterr.perr[i][1];
		if (opt->ionoopt == IONOOPT_IFLC)
		{
			a2 = opt->exterr.perr[i][2];
			b2 = opt->exterr.perr[i][3];
			a = sqrt(SQR(2.55)*a*a + SQR(1.55)*a2*a2);
			b = sqrt(SQR(2.55)*b*b + SQR(1.55)*b2*b2);
		}
	}
	else
	{ /* normal error model */
		if (type == 1) fact *= opt->eratio[0];
		fact *= sys == SYS_GLO ? EFACT_GLO : (sys == SYS_SBS ? EFACT_SBS : EFACT_GPS);
		if (opt->ionoopt == IONOOPT_IFLC) fact *= 3.0;
		a = fact*opt->err[1];
		b = fact*opt->err[2];
	}
	return a*a + b*b / sinel / sinel;
}
/* initialize state and covariance -------------------------------------------*/
static void initx(rtk_t *rtk, double xi, double var, int i)
{
	int j = 0;
	rtk->x[i] = xi;
	for (j = 0; j < rtk->nx; j++)
	{
		if (i == j)
		{
			rtk->P[i + j*rtk->nx] = rtk->P[j + i*rtk->nx] = var;
		}
		else
		{
			rtk->P[i + j*rtk->nx] = rtk->P[j + i*rtk->nx] = 0.0;
		}
	}
}
/* dual-frequency iono-free measurements -------------------------------------*/
static int ifmeas(const obsd_t *obs, const nav_t *nav, const double *azel,
	const prcopt_t *opt, const double *dantr, const double *dants,
	double phw, double *meas, double *var)
{
	const double *lam = nav->lam[obs->sat - 1];
	double c1, c2, L1, L2, P1, P2, P1_C1, P2_C2, gamma;
	int i = 0, j = 1, k;

	trace(4, "ifmeas  :\n");

	/* L1-L2 for GPS/GLO/QZS, L1-L5 for GAL/SBS */
	if (NFREQ >= 3 && (satsys(obs->sat, NULL)&(SYS_GAL | SYS_SBS))) j = 2;

	if (NFREQ < 2 || lam[i] == 0.0 || lam[j] == 0.0) return 0;

	/* test snr mask */
	if (testsnr(0, i, azel[1], obs->SNR[i] * 0.25, &opt->snrmask) ||
		testsnr(0, j, azel[1], obs->SNR[j] * 0.25, &opt->snrmask)) {
		return 0;
	}
	gamma = SQR(lam[j]) / SQR(lam[i]); /* f1^2/f2^2 */
	c1 = gamma / (gamma - 1.0);  /*  f1^2/(f1^2-f2^2) */
	c2 = -1.0 / (gamma - 1.0);  /* -f2^2/(f1^2-f2^2) */

	L1 = obs->L[i] * lam[i]; /* cycle -> m */
	L2 = obs->L[j] * lam[j];
	P1 = obs->P[i];
	P2 = obs->P[j];
	P1_C1 = nav->cbias[obs->sat - 1][1];
	P2_C2 = nav->cbias[obs->sat - 1][2];
	if (opt->sateph == EPHOPT_LEX) {
		P1_C1 = nav->lexeph[obs->sat - 1].isc[0] * CLIGHT; /* ISC_L1C/A */
	}
	if (L1 == 0.0 || L2 == 0.0 || P1 == 0.0 || P2 == 0.0) return 0;

	/* iono-free phase with windup correction */
	meas[0] = c1*L1 + c2*L2 - (c1*lam[i] + c2*lam[j])*phw;

	/* iono-free code with dcb correction */
	if (obs->code[i] == CODE_L1C) P1 += P1_C1; /* C1->P1 */
	if (obs->code[j] == CODE_L2C) P2 += P2_C2; /* C2->P2 */
	meas[1] = c1*P1 + c2*P2;
	var[1] = SQR(ERR_CBIAS);

	if (opt->sateph == EPHOPT_SBAS) meas[1] -= P1_C1; /* sbas clock based C1 */

	/* gps-glonass h/w bias correction for code */
	if (opt->exterr.ena[3] && satsys(obs->sat, NULL) == SYS_GLO) {
		meas[1] += c1*opt->exterr.gpsglob[0] + c2*opt->exterr.gpsglob[1];
	}
	/* antenna phase center variation correction */
	for (k = 0; k < 2; k++) {
		if (dants) meas[k] -= c1*dants[i] + c2*dants[j];
		if (dantr) meas[k] -= c1*dantr[i] + c2*dantr[j];
	}
	return 1;
}
/* get tgd parameter (m) -----------------------------------------------------*/
static double gettgd(int sat, const nav_t *nav)
{
	int i;
	for (i = 0; i < nav->n; i++) {
		if (nav->eph[i].sat != sat) continue;
		return CLIGHT*nav->eph[i].tgd[0];
	}
	return 0.0;
}
/* slant ionospheric delay ---------------------------------------------------*/
static int corr_ion(gtime_t time, const nav_t *nav, int sat, const double *pos,
	const double *azel, int ionoopt, double *ion, double *var,
	int *brk)
{
#ifdef EXTSTEC
	double rate;
#endif
	/* sbas ionosphere model */
	if (ionoopt == IONOOPT_SBAS) {
		return sbsioncorr(time, nav, pos, azel, ion, var);
	}
	/* ionex tec model */
	if (ionoopt == IONOOPT_TEC) {
		return iontec(time, nav, pos, azel, 1, ion, var);
	}
#ifdef EXTSTEC
	/* slant tec model */
	if (ionoopt==IONOOPT_STEC) {
		return stec_ion(time,nav,sat,pos,azel,ion,&rate,var,brk);
	}
#endif
	/* broadcast model */
	if (ionoopt == IONOOPT_BRDC) {
		*ion = ionmodel(time, nav->ion_gps, pos, azel);
		*var = SQR(*ion*ERR_BRDCI);
		return 1;
	}
	/* ionosphere model off */
	*ion = 0.0;
	*var = VAR_IONO_OFF;
	return 1;
}
/* ionosphere and antenna corrected measurements -----------------------------*/
static int corrmeas(const obsd_t *obs, const nav_t *nav, const double *pos,
	const double *azel, const prcopt_t *opt,
	const double *dantr, const double *dants, double phw,
	double *meas, double *var, int *brk)
{
	const double *lam = nav->lam[obs->sat - 1];
	double ion = 0.0, L1, P1, PC, P1_P2, P1_C1, vari, gamma;
	int i;

	trace(4, "corrmeas:\n");

	meas[0] = meas[1] = var[0] = var[1] = 0.0;

	/* iono-free LC */
	if (opt->ionoopt == IONOOPT_IFLC)
	{
		return ifmeas(obs, nav, azel, opt, dantr, dants, phw, meas, var);
	}
	if (lam[0] == 0.0 || obs->L[0] == 0.0 || obs->P[0] == 0.0) return 0;

	if (testsnr(0, 0, azel[1], obs->SNR[0] * 0.25, &opt->snrmask)) return 0;

	L1 = obs->L[0] * lam[0];
	P1 = obs->P[0];

	/* dcb correction */
	gamma = SQR(lam[1] / lam[0]); /* f1^2/f2^2 */
	P1_P2 = nav->cbias[obs->sat - 1][0];
	P1_C1 = nav->cbias[obs->sat - 1][1];
	if (P1_P2 == 0.0 && (satsys(obs->sat, NULL)&(SYS_GPS | SYS_GAL | SYS_QZS))) {
		P1_P2 = (1.0 - gamma)*gettgd(obs->sat, nav);
	}
	if (obs->code[0] == CODE_L1C) P1 += P1_C1; /* C1->P1 */
	PC = P1 - P1_P2 / (1.0 - gamma);               /* P1->PC */

	/* slant ionospheric delay L1 (m) */
	if (!corr_ion(obs->time, nav, obs->sat, pos, azel, opt->ionoopt, &ion, &vari, brk)) {

		trace(2, "iono correction error: time=%s sat=%2d ionoopt=%d\n",
			time_str(obs->time, 2), obs->sat, opt->ionoopt);
		return 0;
	}
	/* ionosphere and windup corrected phase and code */
	meas[0] = L1 + ion - lam[0] * phw;
	meas[1] = PC - ion;

	var[0] += vari;
	var[1] += vari + SQR(ERR_CBIAS);

	/* antenna phase center variation correction */
	for (i = 0; i < 2; i++) {
		if (dants) meas[i] -= dants[0];
		if (dantr) meas[i] -= dantr[0];
	}
	return 1;
}
/* L1/L2 geometry-free phase measurement -------------------------------------*/
static double gfmeas(const obsd_t *obs, const nav_t *nav)
{
	const double *lam = nav->lam[obs->sat - 1];

	if (lam[0] == 0.0 || lam[1] == 0.0 || obs->L[0] == 0.0 || obs->L[1] == 0.0) return 0.0;

	return lam[0] * obs->L[0] - lam[1] * obs->L[1];
}
/* temporal update of position -----------------------------------------------*/
static void udpos_ppp(rtk_t *rtk)
{
	int i;
	trace(3, "udpos_ppp:\n");

	/* fixed mode */
	if (rtk->opt.mode == PMODE_PPP_FIXED)
	{
		for (i = 0; i < 3; i++) initx(rtk, rtk->opt.ru[i], 1E-8, i);
		return;
	}
	/* initialize position for first epoch */
	if (norm(rtk->x, 3) <= 0.0)
	{
		for (i = 0; i < 3; i++) initx(rtk, rtk->sol.rr[i], VAR_POS, i);
	}
	/* static ppp mode */
	if (rtk->opt.mode == PMODE_PPP_STATIC) return;

	/* kinmatic mode without dynamics */
	for (i = 0; i < 3; i++)
	{
		initx(rtk, rtk->sol.rr[i], VAR_POS, i);
	}
}
/* temporal update of clock --------------------------------------------------*/
static void udclk_ppp(rtk_t *rtk)
{
	double dtr;
	int i;

	trace(3, "udclk_ppp:\n");

	/* initialize every epoch for clock (white noise) */
	for (i = 0; i < NSYS; i++)
	{
		if (rtk->opt.sateph == EPHOPT_PREC)
		{
			/* time of prec ephemeris is based gpst */
			/* negelect receiver inter-system bias  */
			dtr = rtk->sol.dtr[0];
		}
		else
		{
			dtr = i == 0 ? rtk->sol.dtr[0] : rtk->sol.dtr[0] + rtk->sol.dtr[i];
		}

		initx(rtk, CLIGHT*dtr, VAR_CLK, IC(i, &rtk->opt));
	}
}
/* temporal update of tropospheric parameters --------------------------------*/
static void udtrop_ppp(rtk_t *rtk)
{
	double pos[3], azel[] = { 0.0, PI / 2.0 }, ztd, var;

	int i = IT(&rtk->opt), j;

	trace(3, "udtrop_ppp:\n");

	if (rtk->x[i] == 0.0)  //第一个历元 
	{
		ecef2pos(rtk->sol.rr, pos);
		ztd = sbstropcorr(rtk->sol.time, pos, azel, &var);
		initx(rtk, ztd, var, i);

		if (rtk->opt.tropopt >= TROPOPT_ESTG)  //估计梯度 
		{
			for (j = 0; j < 2; j++) initx(rtk, 1E-6, VAR_GRA, ++i);
		}
	}
	else
	{
		rtk->P[i*(1 + rtk->nx)] += SQR(rtk->opt.prn[2])*fabs(rtk->tt);

		if (rtk->opt.tropopt >= TROPOPT_ESTG)
		{
			for (j = 0; j < 2; j++)
			{
				rtk->P[++i*(1 + rtk->nx)] += SQR(rtk->opt.prn[2] * 0.1)*fabs(rtk->tt);
			}
		}
	}
}
/* detect cycle slip by LLI --------------------------------------------------*/
static void detslp_ll(rtk_t *rtk, const obsd_t *obs, int n)
{
	int i, j;

	trace(3, "detslp_ll: n=%d\n", n);

	for (i = 0; i < n&&i < MAXOBS; i++) for (j = 0; j < rtk->opt.nf; j++) {
		if (obs[i].L[j] == 0.0 || !(obs[i].LLI[j] & 3)) continue;

		trace(3, "detslp_ll: slip detected sat=%2d f=%d\n", obs[i].sat, j + 1);

		rtk->ssat[obs[i].sat - 1].slip[j] = 1;
	}
}
/* detect cycle slip by geometry free phase jump -----------------------------*/
static void detslp_gf(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
	double g0, g1;
	int i, j;

	trace(3, "detslp_gf: n=%d\n", n);

	for (i = 0; i < n&&i<MAXOBS; i++)
	{

		if ((g1 = gfmeas(obs + i, nav)) == 0.0) continue;

		g0 = rtk->ssat[obs[i].sat - 1].gf;
		rtk->ssat[obs[i].sat - 1].gf = g1;

		trace(4, "detslip_gf: sat=%2d gf0=%8.3f gf1=%8.3f\n", obs[i].sat, g0, g1);

		if (g0 != 0.0&&fabs(g1 - g0)>rtk->opt.thresslip)
		{
			trace(3, "detslip_gf: slip detected sat=%2d gf=%8.3f->%8.3f\n",
				obs[i].sat, g0, g1);

			for (j = 0; j < rtk->opt.nf; j++) rtk->ssat[obs[i].sat - 1].slip[j] |= 1;
		}
	}
}

/* temporal update of phase biases -------------------------------------------*/
static void udbias_ppp(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
	double meas[2], var[2], bias[MAXOBS] = { 0 }, offset = 0.0, pos[3] = { 0 };
	int i, j, k, sat, brk = 0;

	trace(3, "udbias  : n=%d\n", n);

	for (i = 0; i < MAXSAT; i++) for (j = 0; j < rtk->opt.nf; j++)
	{
		rtk->ssat[i].slip[j] = 0;
	}
	/* detect cycle slip by LLI */
	detslp_ll(rtk, obs, n);

	/* detect cycle slip by geometry-free phase jump */
	detslp_gf(rtk, obs, n, nav);

	/* reset phase-bias if expire obs outage counter */
	for (i = 0; i<MAXSAT; i++)
	{
		if (++rtk->ssat[i].outc[0]>(unsigned int)rtk->opt.maxout)
		{
			initx(rtk, 0.0, 0.0, IB(i + 1, &rtk->opt));
		}
	}

	ecef2pos(rtk->sol.rr, pos);

	for (i = k = 0; i < n&&i < MAXOBS; i++)
	{
		sat = obs[i].sat;
		j = IB(sat, &rtk->opt);
		if (!corrmeas(obs + i, nav, pos, rtk->ssat[sat - 1].azel, &rtk->opt, NULL, NULL,
			0.0, meas, var, &brk)) continue;

		if (brk)
		{
			rtk->ssat[sat - 1].slip[0] = 1;
			trace(2, "%s: sat=%2d correction break\n", time_str(obs[i].time, 0), sat);
		}

		bias[i] = meas[0] - meas[1];

		if (rtk->x[j] == 0.0 ||
			rtk->ssat[sat - 1].slip[0] || rtk->ssat[sat - 1].slip[1]) continue;
		offset += bias[i] - rtk->x[j];
		k++;
	}
	/* correct phase-code jump to enssure phase-code coherency */
	if (k >= 2 && fabs(offset / k) > 0.0005*CLIGHT)
	{
		for (i = 0; i < MAXSAT; i++)
		{
			j = IB(i + 1, &rtk->opt);
			if (rtk->x[j] != 0.0) rtk->x[j] += offset / k;
		}

		trace(2, "phase-code jump corrected: %s n=%2d dt=%12.9fs\n",
			time_str(rtk->sol.time, 0), k, offset / k / CLIGHT);
	}
	for (i = 0; i < n&&i < MAXOBS; i++)
	{
		sat = obs[i].sat;
		j = IB(sat, &rtk->opt);

		rtk->P[j + j*rtk->nx] += SQR(rtk->opt.prn[0])*fabs(rtk->tt);

		if (rtk->x[j] != 0.0&&
			!rtk->ssat[sat - 1].slip[0] && !rtk->ssat[sat - 1].slip[1]) continue;

		if (bias[i] == 0.0) continue;

		/* reinitialize phase-bias if detecting cycle slip */
		initx(rtk, bias[i], VAR_BIAS, IB(sat, &rtk->opt));

		trace(5, "udbias_ppp: sat=%2d bias=%.3f\n", sat, meas[0] - meas[1]);
	}
}



/* temporal update of states --------------------------------------------------*/
static void udstate_ppp(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
	trace(3, "udstate_ppp: n=%d\n", n);

	/* temporal update of position */
	udpos_ppp(rtk);  //前3个是XYZ坐标

	/* temporal update of clock */
	udclk_ppp(rtk);  //每个卫星系统估计相对于GPS系统的卫星钟差，既接收机钟差+系统时延

	/* temporal update of tropospheric parameters */
	if (rtk->opt.tropopt >= TROPOPT_EST)
	{
		udtrop_ppp(rtk);
	}
	/* temporal update of phase-bias */
	udbias_ppp(rtk, obs, n, nav);
}
/* satellite antenna phase center variation ----------------------------------*/
static void satantpcv(const double *rs, const double *rr, const pcv_t *pcv,
	double *dant)
{
	double ru[3], rz[3], eu[3], ez[3], nadir, cosa;
	int i;

	for (i = 0; i < 3; i++) {
		ru[i] = rr[i] - rs[i];
		rz[i] = -rs[i];
	}
	if (!normv3(ru, eu) || !normv3(rz, ez)) return;

	cosa = dot(eu, ez, 3);
	cosa = cosa<-1.0 ? -1.0 : (cosa>1.0 ? 1.0 : cosa);
	nadir = acos(cosa);

	antmodel_s(pcv, nadir, dant);
}
/* precise tropospheric model ------------------------------------------------*/
static double prectrop(gtime_t time, const double *pos, const double *azel,
	const prcopt_t *opt, const double *x, double *dtdx,
	double *var)
{
	const double zazel[] = { 0.0, PI / 2.0 };
	double zhd, m_h, m_w, cotz, grad_n, grad_e;

	/* zenith hydrostatic delay */
	zhd = tropmodel(time, pos, zazel, 0.0);

	/* mapping function */
	m_h = tropmapf(time, pos, azel, &m_w);

	if ((opt->tropopt == TROPOPT_ESTG || opt->tropopt == TROPOPT_CORG) && azel[1] > 0.0) {

		/* m_w=m_0+m_0*cot(el)*(Gn*cos(az)+Ge*sin(az)): ref [6] */
		cotz = 1.0 / tan(azel[1]);
		grad_n = m_w*cotz*cos(azel[0]);
		grad_e = m_w*cotz*sin(azel[0]);
		m_w += grad_n*x[1] + grad_e*x[2];
		dtdx[1] = grad_n*(x[0] - zhd);
		dtdx[2] = grad_e*(x[0] - zhd);
	}
	dtdx[0] = m_w;
	*var = SQR(0.01);
	return m_h*zhd + m_w*(x[0] - zhd);
}
/* phase and code residuals --------------------------------------------------*/
static int res_ppp(int iter, const obsd_t *obs, int n, const double *rs,
	const double *dts, const double *vare, const int *svh,
	const nav_t *nav, const double *x, rtk_t *rtk, double *v,
	double *H, double *R, double *azel)
{
	prcopt_t *opt = &rtk->opt;
	double r, rr[3], disp[3], pos[3], e[3], meas[2], dtdx[3], dantr[NFREQ] = { 0 };
	double dants[NFREQ] = { 0 }, var[MAXOBS * 2], dtrp = 0.0, vart = 0.0, varm[2] = { 0 };
	int i, j, k, sat, sys, nv = 0, nx = rtk->nx, brk, tideopt;

	trace(3, "res_ppp : n=%d nx=%d\n", n, nx);

	for (i = 0; i < MAXSAT; i++) rtk->ssat[i].vsat[0] = 0;

	for (i = 0; i < 3; i++) rr[i] = x[i];

	/* earth tides correction */
	if (opt->tidecorr) {
		tideopt = opt->tidecorr == 1 ? 1 : 7; /* 1:solid, 2:solid+otl+pole */

		tidedisp(gpst2utc(obs[0].time), rr, tideopt, &nav->erp, opt->odisp[0],
			disp);
		for (i = 0; i < 3; i++) rr[i] += disp[i];
	}
	ecef2pos(rr, pos);

	for (i = 0; i < n&&i < MAXOBS; i++)
	{
		sat = obs[i].sat;

		if (!(sys = satsys(sat, NULL)) || !rtk->ssat[sat - 1].vs) continue;

		/* geometric distance/azimuth/elevation angle */
		if ((r = geodist(rs + i * 6, rr, e)) <= 0.0 ||
			satazel(pos, e, azel + i * 2) < opt->elmin) continue;

		/* excluded satellite? */
		if (satexclude(obs[i].sat, svh[i], opt)) continue;

		/* tropospheric delay correction */
		if (opt->tropopt == TROPOPT_SAAS)
		{
			dtrp = tropmodel(obs[i].time, pos, azel + i * 2, REL_HUMI);
			vart = SQR(ERR_SAAS);
		}
		else if (opt->tropopt == TROPOPT_SBAS)
		{
			dtrp = sbstropcorr(obs[i].time, pos, azel + i * 2, &vart);
		}
		else if (opt->tropopt == TROPOPT_EST || opt->tropopt == TROPOPT_ESTG)
		{
			dtrp = prectrop(obs[i].time, pos, azel + i * 2, opt, x + IT(opt), dtdx, &vart);
		}
		else if (opt->tropopt == TROPOPT_COR || opt->tropopt == TROPOPT_CORG)
		{
			dtrp = prectrop(obs[i].time, pos, azel + i * 2, opt, x, dtdx, &vart);
		}
		/* satellite antenna model */
		if (opt->posopt[0])
		{
			satantpcv(rs + i * 6, rr, nav->pcvs + sat - 1, dants);
		}
		/* receiver antenna model */
		antmodel(opt->pcvr, opt->antdel[0], azel + i * 2, opt->posopt[1], dantr);

		/* phase windup correction */
		if (opt->posopt[2])
		{
			windupcorr(rtk->sol.time, rs + i * 6, rr, &rtk->ssat[sat - 1].phw);
		}
		/* ionosphere and antenna phase corrected measurements */
		if (!corrmeas(obs + i, nav, pos, azel + i * 2, &rtk->opt, dantr, dants,
			rtk->ssat[sat - 1].phw, meas, varm, &brk))
		{
			continue;
		}
		/* satellite clock and tropospheric delay */
		r += -CLIGHT*dts[i * 2] + dtrp;

		trace(5, "sat=%2d azel=%6.1f %5.1f dtrp=%.3f dantr=%6.3f %6.3f dants=%6.3f %6.3f phw=%6.3f\n",
			sat, azel[i * 2] * R2D, azel[1 + i * 2] * R2D, dtrp, dantr[0], dantr[1], dants[0],
			dants[1], rtk->ssat[sat - 1].phw);

		for (j = 0; j < 2; j++)
		{ /* for phase and code */

			if (meas[j] == 0.0) continue;

			for (k = 0; k < nx; k++) H[k + nx*nv] = 0.0;

			v[nv] = meas[j] - r;

			for (k = 0; k < 3; k++) H[k + nx*nv] = -e[k];

			if (sys != SYS_GLO)
			{
				v[nv] -= x[IC(0, opt)];
				H[IC(0, opt) + nx*nv] = 1.0;
			}
			else
			{
				v[nv] -= x[IC(1, opt)];
				H[IC(1, opt) + nx*nv] = 1.0;
			}
			if (opt->tropopt >= TROPOPT_EST)
			{
				for (k = 0; k < (opt->tropopt >= TROPOPT_ESTG ? 3 : 1); k++)
				{
					H[IT(opt) + k + nx*nv] = dtdx[k];
				}
			}
			if (j == 0)
			{
				v[nv] -= x[IB(obs[i].sat, opt)];
				H[IB(obs[i].sat, opt) + nx*nv] = 1.0;
			}
			var[nv] = varerr(obs[i].sat, sys, azel[1 + i * 2], j, opt) + varm[j] + vare[i] + vart;

			if (j == 0) rtk->ssat[sat - 1].resc[0] = v[nv];
			else      rtk->ssat[sat - 1].resp[0] = v[nv];

			/* test innovation */
#if 0
			if (opt->maxinno>0.0&&fabs(v[nv])>opt->maxinno) 
			{
#else
			if (opt->maxinno > 0.0&&fabs(v[nv]) > opt->maxinno&&sys != SYS_GLO)
			{
#endif
				trace(2, "ppp outlier rejected %s sat=%2d type=%d v=%.3f\n",
					time_str(obs[i].time, 0), sat, j, v[nv]);
				rtk->ssat[sat - 1].rejc[0]++;
				continue;
			}
			if (j == 0) rtk->ssat[sat - 1].vsat[0] = 1;
			nv++;
		}
	}
	for (i = 0; i < nv; i++) for (j = 0; j < nv; j++) {
		R[i + j*nv] = i == j ? var[i] : 0.0;
	}
	trace(5, "x=\n"); tracemat(5, x, 1, nx, 8, 3);
	trace(5, "v=\n"); tracemat(5, v, 1, nv, 8, 3);
	trace(5, "H=\n"); tracemat(5, H, nx, nv, 8, 3);
	trace(5, "R=\n"); tracemat(5, R, nv, nv, 8, 5);
	return nv;
}
/* number of estimated states ------------------------------------------------*/
extern int pppnx(const prcopt_t *opt)
{
	return NX(opt);
}
/* precise point positioning -------------------------------------------------*/
extern void pppos(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
	const prcopt_t *opt = &rtk->opt;
	double *rs, *dts, *var, *v, *H, *R, *azel, *xp, *Pp;
	int i, nv, info, svh[MAXOBS], stat = SOLQ_SINGLE;

	trace(3, "pppos   : nx=%d n=%d\n", rtk->nx, n);

	rs = mat(6, n); dts = mat(2, n); var = mat(1, n); azel = zeros(2, n);

	for (i = 0; i < MAXSAT; i++) rtk->ssat[i].fix[0] = 0;

	/* temporal update of states */
	udstate_ppp(rtk, obs, n, nav);

	trace(4, "x(0)="); tracemat(4, rtk->x, 1, NR(opt), 13, 4);

	/* satellite positions and clocks */
	satposs(obs[0].time, obs, n, nav, rtk->opt.sateph, rs, dts, var, svh);

	/* exclude measurements of eclipsing satellite */
	if (rtk->opt.posopt[3])
	{
		testeclipse(obs, n, nav, rs);
	}
	xp = mat(rtk->nx, 1); Pp = zeros(rtk->nx, rtk->nx);
	matcpy(xp, rtk->x, rtk->nx, 1);
	nv = n*rtk->opt.nf * 2; v = mat(nv, 1); H = mat(rtk->nx, nv); R = mat(nv, nv);

	for (i = 0; i < rtk->opt.niter; i++)
	{

		/* phase and code residuals */
		if ((nv = res_ppp(i, obs, n, rs, dts, var, svh, nav, xp, rtk, v, H, R, azel)) <= 0) break;

		/* measurement update */
		matcpy(Pp, rtk->P, rtk->nx, rtk->nx);

		if ((info = filter(xp, Pp, H, v, R, rtk->nx, nv)))
		{
			trace(2, "ppp filter error %s info=%d\n", time_str(rtk->sol.time, 0),
				info);
			break;
		}
		trace(4, "x(%d)=", i + 1); tracemat(4, xp, 1, NR(opt), 13, 4);

		stat = SOLQ_PPP;
	}
	if (stat == SOLQ_PPP)
	{
		/* postfit residuals */
		res_ppp(1, obs, n, rs, dts, var, svh, nav, xp, rtk, v, H, R, azel);

		/* update state and covariance matrix */
		matcpy(rtk->x, xp, rtk->nx, 1);
		matcpy(rtk->P, Pp, rtk->nx, rtk->nx);

		/* ambiguity resolution in ppp */
		if (opt->modear == ARMODE_PPPAR || opt->modear == ARMODE_PPPAR_ILS) {
			if (pppamb(rtk, obs, n, nav, azel)) stat = SOLQ_FIX;
		}
		/* update solution status */
		rtk->sol.ns = 0;
		for (i = 0; i < n&&i < MAXOBS; i++) {
			if (!rtk->ssat[obs[i].sat - 1].vsat[0]) continue;
			rtk->ssat[obs[i].sat - 1].lock[0]++;
			rtk->ssat[obs[i].sat - 1].outc[0] = 0;
			rtk->ssat[obs[i].sat - 1].fix[0] = 4;
			rtk->sol.ns++;
		}
		rtk->sol.stat = stat;

		for (i = 0; i < 3; i++) {
			rtk->sol.rr[i] = rtk->x[i];
			rtk->sol.qr[i] = (float)rtk->P[i + i*rtk->nx];
		}
		rtk->sol.qr[3] = (float)rtk->P[1];
		rtk->sol.qr[4] = (float)rtk->P[2 + rtk->nx];
		rtk->sol.qr[5] = (float)rtk->P[2];
		rtk->sol.dtr[0] = rtk->x[IC(0, opt)];
		rtk->sol.dtr[1] = rtk->x[IC(1, opt)] - rtk->x[IC(0, opt)];
		for (i = 0; i < n&&i < MAXOBS; i++) {
			rtk->ssat[obs[i].sat - 1].snr[0] = MIN(obs[i].SNR[0], obs[i].SNR[1]);
		}
		for (i = 0; i < MAXSAT; i++) {
			if (rtk->ssat[i].slip[0] & 3) rtk->ssat[i].slipc[0]++;
		}
	}
	free(rs); free(dts); free(var); free(azel);
	free(xp); free(Pp); free(v); free(H); free(R);
}


/*自定义的状态更新函数*/
extern int myUdstate_ppp(rtk_t* rtk, const obsd_t* obs, int n, const nav_t* nav)
{
	int stat = 0, stat1 = 0, stat2 = 0, j = 0;
	trace(3, "udstate_ppp: n=%d\n", n);

	/* temporal update of position */
	stat1 = myUdpos_ppp(rtk);  //前3个是XYZ坐标

	/* temporal update of clock */
	stat2 = myUdclk_ppp(rtk);  //每个卫星系统估计相对于GPS系统的卫星钟差，既接收机钟差+系统时延

	/* temporal update of tropospheric parameters */
	if (rtk->opt.tropopt >= TROPOPT_EST)
	{
		//static void myudtrop_ppp( rtk_t *rtk )
		stat = myUdtrop_ppp(rtk);
	}

	/* temporal update of phase-bias */
	myUpbias_ppp(rtk, obs, n, nav);
}


#define myIB(s,f,opt) ( NR(opt)+(s) - 1 + (f*MAXSAT) )    /* state index of phase bias, s是卫是几号卫星（从1计数），f是第几个频率(从0计数) */

/* initialize rtk control ------------------------------------------------------
* initialize rtk control struct
* args   : rtk_t    *rtk    IO  rtk control/result struct
*          prcopt_t *opt    I   positioning options (see rtklib.h)
* return : none
*-----------------------------------------------------------------------------
Author: Zhen.Li
*/
extern void myrtkinit(rtk_t *rtk, const prcopt_t *opt)
{
	sol_t sol0 = { { 0 } };
	ambc_t ambc0 = { { { 0 } } };
	ssat_t ssat0 = { 0 };
	int i;

	trace(3, "rtkinit :\n");

	rtk->sol = sol0;
	for (i = 0; i < 6; i++) rtk->rb[i] = 0.0;
	//rtk->nx=opt->mode<=PMODE_FIXED?NX(opt):pppnx(opt);
	rtk->nx = NR(opt) + MAXSAT*opt->nf;
	rtk->na = NR(opt);  //除去模糊度以外，其他的参数个数

	rtk->tt = 0.0;
	rtk->x = zeros(rtk->nx, 1);
	rtk->P = zeros(rtk->nx, rtk->nx);
	rtk->xa = zeros(rtk->na, 1);
	rtk->Pa = zeros(rtk->na, rtk->na);
	rtk->nfix = rtk->neb = 0;
	for (i = 0; i < MAXSAT; i++)
	{
		rtk->ambc[i] = ambc0;
		rtk->ssat[i] = ssat0;
	}
	for (i = 0; i < MAXERRMSG; i++) rtk->errbuf[i] = 0;
	rtk->opt = *opt;
}


/*自定义的模糊度和周跳处理函数
Author: Zhen.Li
*/
extern int myUpbias_ppp(rtk_t* rtk, const obsd_t* obs, int satnum, const nav_t* nav)
{

	double* myx, *myp;    //注意释放内存


	//double meas[2];
	double var[2] = { 0.0 }, bias[MAXSAT*NFREQ] = { 0 }, offset = 0.0, pos[3] = { 0 }, ion = 0;
	int i, j, k, s, sat, brk = 0, f = -1, indexC = 0, indexT = 0;
	double test1 = 0.0, test2 = 0.0, test = 0.0;

	trace(3, "udbias  : n=%d\n", satnum);


	myx = mat(rtk->nx, 1);   //rtk->nx是所有待估计的参数个数 
	myp = mat(rtk->nx, rtk->nx);
	memset(myx, 0, sizeof(double)*rtk->nx);
	memset(myp, 0, sizeof(double)*rtk->nx*rtk->nx);


	matcpy(myx, rtk->x, rtk->nx, 1); // rtk->x 是状态向量
	matcpy(myp, rtk->P, rtk->nx, rtk->nx);

	memset(rtk->x, 0, sizeof(double)*rtk->nx);  //对rtk->x进行清0
	//memset(rtk->P,0,sizeof(double)*rtk->nx*rtk->nx);// 对rtk->p进行清0

	//拷贝坐标
	for (i = 0; i < 3; i++)
	{
		//initx(rtk,myx[i],myp[i*rtk->nx+i],i);
		rtk->x[i] = myx[i];
	}

	indexC = IC(0, &rtk->opt);
	indexT = IT(&rtk->opt);
	rtk->x[indexC] = myx[indexC];
	rtk->x[indexT] = myx[indexT];

	//initx(rtk,myx[indexC],myp[indexC*(rtk->nx)+indexC],indexC);
	//initx(rtk,myx[indexT],myp[indexT*(rtk->nx)+indexT],indexT );


	fprintf(g_mylog, "before update_bias\n");
	fprintf(g_mylog, "x :");
	for (i = k = 0; i < satnum&&i < MAXOBS; i++)  //这里到底是MAXOBS还是MAXSAT呢? 
	{
		sat = obs[i].sat;
		for (f = 0; f < rtk->opt.nf; f++)  //对所有的频率循环
		{
			j = myIB(sat, f, &rtk->opt);  //当前卫星，当前频率的模糊度的序号
			if (rtk->x[j] != 0.0)
			{
				fprintf(g_mylog, "%5.3f  ", rtk->x[j]);
			}
		}
	}
	fprintf(g_mylog, "\n");


	//周跳标志清0
	for (i = 0; i < MAXSAT; i++)
	{
		for (j = 0; j < rtk->opt.nf; j++)
		{
			rtk->ssat[i].slip[j] = 0;
		}
	}

	/* detect cycle slip by LLI */
	detslp_ll(rtk, obs, satnum);  //这个从观测值中获取当前历元是否发生了周跳

	/* detect cycle slip by geometry-free phase jump */
	detslp_gf(rtk, obs, satnum, nav);  //双频只探测不修复

	/*进行周跳探测*/



	///* reset phase-bias if expire obs outage counter */
	//for ( i=0;i<MAXSAT;i++ )
	//{
	//		for( f =0; f< rtk->opt.nf ; f++ )
	//		{
	//			if ( ++rtk->ssat[i].outc[f]>(unsigned int)rtk->opt.maxout)  //这里只用考虑第1个频率
	//			{
	//				k = myIB(i+1,f,&rtk->opt);
	//				initx(rtk,0.0,0.0,myIB(i+1,f,&rtk->opt));
	//			}	
	//		}	
	//}

	ecef2pos(rtk->sol.rr, pos);

	for (i = k = 0; i < satnum&&i < MAXOBS; i++)  //这里到底是MAXOBS还是MAXSAT呢? 
	{
		sat = obs[i].sat;
		for (f = 0; f < rtk->opt.nf; f++)  //对所有的频率循环
		{
			j = myIB(sat, f, &rtk->opt);  //当前卫星，当前频率的模糊度的序号

			test1 = (obs + i)->P[f];   //伪距											
			test2 = (obs + i)->L[f];   //载波

			if (fabs(test1) <= 100.0 || fabs(test2) <= 100.0)  continue;

			corr_ion(obs->time, nav, sat, pos, rtk->ssat[sat - 1].azel, rtk->opt.ionoopt, &ion, var, brk);

			//设置模糊初始值(应该是整周)
			bias[k*satnum + f] = ((test2*nav->lam[sat - 1][f]) - (test1)) / nav->lam[sat - 1][f];  //模糊的初始值，单位：周

			if (myx[j] == 0.0
				|| rtk->ssat[sat - 1].slip[0]
				|| rtk->ssat[sat - 1].slip[1]
				)
			{
				continue;
			}

			offset += (bias[k*satnum + f] - myx[j]);  //记录下特大周跳 unit:cycle
		}

		k++;
	}


	/* correct phase-code jump to enssure phase-code coherency */
	//if (k>=2&&fabs(offset/k)>0.0005*CLIGHT*0.2 ) 
	//{
	//	for ( i=0;i<MAXSAT;i++)
	//	{
	//		for(f =0; f<rtk->opt.nf; f++ )
	//		{
	//			j = myIB(i+1,f,&rtk->opt);
	//			if ( fabs(bias[i*satnum+f])> 500 )  //500周以上的算大周跳 
	//			{
	//				rtk->x[j]=  bias[i*satnum+f];  //调整初始模糊度的值
	//				rtk->P[j+j*rtk->nx] = 100*100;
	//			}
	//		}	
	//	}
	//	trace(2,"phase-code jump corrected: %s n=%2d dt=%12.9fs\n",
	//		time_str(rtk->sol.time,0),k,offset/k/CLIGHT);
	//}



	for (i = k = 0; i < satnum&&i < MAXOBS; i++)  //这里到底是MAXOBS还是MAXSAT呢? 
	{
		sat = obs[i].sat;
		for (f = 0; f < rtk->opt.nf; f++)  //对所有的频率循环
		{
			j = myIB(sat, f, &rtk->opt);  //当前卫星，当前频率的模糊度的序号

			//rtk->P[j+j*rtk->nx] = myp[j+j*rtk->nx] +  SQR(rtk->opt.prn[0])*fabs(rtk->tt);

			rtk->P[j + j*rtk->nx] += SQR(rtk->opt.prn[0])*fabs(rtk->tt);

			if (myx[j] != 0.0&&
				!rtk->ssat[sat - 1].slip[0] && !rtk->ssat[sat - 1].slip[1]
				)
			{
				rtk->x[j] = myx[j];
				continue;  //上一历元有值，而且没有发生周跳，则保持原值
			}

			if (bias[k*satnum + f] == 0.0 && myx[j] == 0.0)   //上下两个历元均不存在这个卫星，则认为当前卫星不存在，则直接下一个 
			{
				rtk->x[j] = myx[j];
				continue;  //第i颗卫星,第j个频率
			}

			//if ( bias[k*satnum+f]==0.0 && rtk->x[j] != 0.0 )   //说明当前历元，这个卫星不存在了
			//{
			//	rtk->x[j] = 0.0;  //清空不存在卫星的信息
			//	initx(rtk,0.0,0.0,myIB(sat,f,&rtk->opt));  //清空不存在卫星的信息
			//	continue;  //第i颗卫星,第j个频率
			//}

			/* reinitialize phase-bias if detecting cycle slip
			//第一个历元也用这个初始化
			*/
			if (rtk->ssat[sat - 1].slip[0] || rtk->ssat[sat - 1].slip[1])
			{
				//initx(rtk,bias[k*n+f],VAR_BIAS,myIB(sat,f,&rtk->opt));
				fprintf(g_mylog, "CYCLE SLIPS : sat -> %d\n", sat);
				//continue;
			}

			initx(rtk, bias[k*satnum + f], VAR_BIAS, myIB(sat, f, &rtk->opt));

			/*myx[j] = bias[k*satnum+f];
			for( s=0;s<rtk->nx;s++ )
			{
			if( s== j )
			{
			myp[j+s*rtk->nx]=rtk->P[s+j*rtk->nx] = VAR_BIAS;
			}
			else
			{
			myp[j+s*rtk->nx]=rtk->P[s+j*rtk->nx] = 0.0;
			}
			}*/



		}

		k++;
	}


	fprintf(g_mylog, "x and bias check:\n");
	fprintf(g_mylog, "x :");
	for (i = k = 0; i < satnum&&i < MAXOBS; i++)  //这里到底是MAXOBS还是MAXSAT呢? 
	{
		sat = obs[i].sat;
		for (f = 0; f < rtk->opt.nf; f++)  //对所有的频率循环
		{
			j = myIB(sat, f, &rtk->opt);  //当前卫星，当前频率的模糊度的序号
			if (rtk->x[j] != 0.0)
			{
				fprintf(g_mylog, "%5.3f  ", rtk->x[j]);
			}
		}
	}
	fprintf(g_mylog, "\n");

	fprintf(g_mylog, "bias :");
	for (i = k = 0; i < satnum&&i < MAXOBS; i++)  //这里到底是MAXOBS还是MAXSAT呢? 
	{
		sat = obs[i].sat;
		for (f = 0; f < rtk->opt.nf; f++)  //对所有的频率循环
		{
			if (bias[k*satnum + f] != 0.0)
			{
				fprintf(g_mylog, "%5.3f  ", bias[k*satnum + f]);
			}
		}
		k++;
	}
	fprintf(g_mylog, "\n");



	fprintf(g_mylog, "after the update :\n");
	for (j = 0; j < rtk->nx; j++)
	{
		if (rtk->x[j] != 0.0)
		{
			fprintf(g_mylog, "%10.3f ", rtk->x[j]);
		}
	}
	fprintf(g_mylog, "\n");


	if (myx != NULL) free(myx);
	if (myp != NULL) free(myp);


}



/* temporal update of position -----------------------------------------------
Author: Zhen.Li
*/
static int myUdpos_ppp(rtk_t *rtk)
{
	int i = 0;
	//trace(3,"udpos_ppp:\n");

	///* fixed mode 固定解初始化*/
	//if( rtk->opt.mode==PMODE_PPP_FIXED )
	//{
	//	for( i = 0;i<3;i++ ) 
	//	{
	//	   initx(rtk,rtk->opt.ru[i],1E-8,i);  //针对于流动站,初始化的对象时rtk->opt.ru
	//	}
	//	return;
	//}
	/* initialize position for first epoch */
	//if ( norm(rtk->x,3)<=0.0 )
	//{
	//	for( i = 0;i< 3; i++ )
	//	{
	//		initx( rtk,rtk->sol.rr[i],VAR_POS,i ); //初始化的对象时rtk->sol.rr
	//	}
	//}
	/* static ppp mode */
	//if ( rtk->opt.mode==PMODE_PPP_STATIC ) return;

	/* kinmatic mode without dynamics 动态模式，每个历元，坐标的方差都重新设置,坐标值保持不变*/
	for (i = 0; i < 3; i++)
	{
		initx(rtk, rtk->sol.rr[i], VAR_POS, i);
	}

	return 0;
}


/* temporal update of clock --------------------------------------------------
Author: Zhen.Li
*/
static int myUdclk_ppp(rtk_t *rtk)
{
	double dtr;

	int index;

	trace(3, "udclk_ppp:\n");

	dtr = rtk->sol.dtr[0];

	dtr = dtr * CLIGHT;

	index = IC(0, &rtk->opt);

	printf("IC_index--%d  %.3f %.3f \n", index, rtk->x[index], rtk->P[index*rtk->nx + index]);

	initx(rtk, dtr, VAR_CLK, index); //VAR_CLK

	/*for (i=0;i<NSYS;i++ )
	{
	dtr=rtk->sol.dtr[0];
	initx(rtk,CLIGHT*dtr,VAR_CLK,IC(i,&rtk->opt));
	}*/

	///* initialize every epoch for clock (white noise) */
	//for (i=0;i<NSYS;i++)
	//{
	//	
	//	dtr=i==0?rtk->sol.dtr[0]:rtk->sol.dtr[0]+rtk->sol.dtr[i];
	//	
	//	initx(rtk,CLIGHT*dtr,VAR_CLK,IC(i,&rtk->opt));
	//}

	return 0;
}

/* temporal update of tropospheric parameters --------------------------------
Author: Zhen.Li
*/
static int myUdtrop_ppp(rtk_t *rtk)
{
	double pos[3], azel[] = { 0.0, PI / 2.0 }, ztd, var, zwd;

	int index = IT(&rtk->opt), j;

	trace(3, "udtrop_ppp:\n");

	if (rtk->x[index] == 0.0)  //第一个历元 
	{
		ecef2pos(rtk->sol.rr, pos);
		zwd = 0.15;
		var = SQR(0.5);
		//ztd=sbstropcorr(rtk->sol.time,pos,azel,&var);

		initx(rtk, zwd, var, index);

		//if (rtk->opt.tropopt>=TROPOPT_ESTG)  //估计梯度 
		//{
		//	for (j=0;j<2;j++) initx(rtk,1E-6,VAR_GRA,++i);
		//}
	}
	else
	{
		//添加过程噪声
		rtk->P[index*(rtk->nx) + index] += SQR(rtk->opt.prn[2])*fabs(rtk->tt);
		//rtk->P[ index*(rtk->nx)+ index ] = 0.5*0.5;

		/*if (rtk->opt.tropopt>=TROPOPT_ESTG)
		{
		for (j=0;j<2;j++)
		{
		rtk->P[++i*(1+rtk->nx)]+=SQR(rtk->opt.prn[2]*0.1)*fabs(rtk->tt);
		}
		}*/
	}
}


/*自定义的PPP非差参数计算
Author: Zhen.Li
*/
extern int myPPP_RES(int iter, const obsd_t *obs, int satnum, const double *rs,
	const double *dts, const double *vare, const int *svh,
	const nav_t *nav, const double *x, rtk_t *rtk, double *v,
	double *H, double *R, double *azel)
{
	const double ep[] = { 2000, 1, 1, 12, 0, 0 };
	prcopt_t *opt = &rtk->opt;
	double geodis, r, rr[3], disp[3], pos[3], e[3], meas[2], dtdx[3], dantr[NFREQ] = { 0 };  //means可能需要调整
	double dants[NFREQ] = { 0 }, var[MAXOBS * 2], dtrp = 0.0, vart = 0.0, varm[2] = { 0 };
	int i, j, k, f, sat, sys, nv = 0, nx = rtk->nx, brk, tideopt;
	double clk_sat = 0.0, clk_sta = 0.0;

	double gamma, P1_P2, P1_C1, alfa = 0.0, beta = 0.0;

	double test1 = 0.0, test2 = 0.0, test = 0.0;

	double zhd, m_h, m_w, cotz, grad_n, grad_e;
	double zazel[2] = { 0.0, PI / 2 };
	int frequenCount = 0;   //对每颗卫星的频率进行计数，保证均为双频载波和双频伪距时才引入载波无电离层观测值

	//当前的mjd计算，方便调用gmf
	double mjd = 51544.5 + (timediff(obs->time, epoch2time(ep))) / 86400.0;

	trace(3, "res_ppp : n=%d nx=%d\n", satnum, nx);

	//注意这里设定了vsat的值为0
	for (i = 0; i < MAXSAT; i++) rtk->ssat[i].vsat[0] = 0;

	for (i = 0; i < 3; i++) rr[i] = x[i];  //复制坐标参数

	/* earth tides correction */
	if (opt->tidecorr)  //是否进行固体潮改正 
	{
		tideopt = opt->tidecorr == 1 ? 1 : 7; /* 1:solid, 2:solid+otl+pole */

		tidedisp(gpst2utc(obs[0].time), rr, tideopt, &nav->erp, opt->odisp[0], disp);
		//添加地球固体潮改正到测站坐标上
		for (i = 0; i < 3; i++) rr[i] += disp[i];
	}

	ecef2pos(rr, pos);

	fprintf(g_mylog, "observable(satnum:%02d):\n", satnum);

	for (i = 0; i < satnum&&i < MAXOBS; i++)  // for all the satellites 
	{
		sat = obs[i].sat;

		if (!(sys = satsys(sat, NULL))  /* || !rtk->ssat[sat-1].vs*/)   //搞清楚这里的vs的具体含义 
		{
			for (f = 0; f < rtk->opt.nf; f++)
			{
				j = myIB(sat, f, &rtk->opt);
				rtk->x[j] = 0.0;  //模糊度设置为0
				rtk->P[j + j*rtk->nx] = 0.0;
			}
			f = 0;
			j = 0;
			continue;
		}

		/* geometric distance/azimuth/elevation angle */
		if ((geodis = geodist(rs + i * 6, rr, e)) <= 0.0 ||
			satazel(pos, e, azel + i * 2) < opt->elmin)
		{
			for (f = 0; f < rtk->opt.nf; f++)
			{
				j = myIB(sat, f, &rtk->opt);
				rtk->x[j] = 0.0;  //模糊度设置为0
				rtk->P[j + j*rtk->nx] = 0.0;
			}
			f = 0; j = 0;
			continue;
		}

		/* excluded satellite? */
		if (satexclude(obs[i].sat, svh[i], opt))
		{
			for (f = 0; f < rtk->opt.nf; f++)
			{
				j = myIB(sat, f, &rtk->opt);
				rtk->x[j] = 0.0;  //模糊度设置为0
				rtk->P[j + j*rtk->nx] = 0.0;
			}
			f = 0;
			j = 0;
			continue;
		}
		/* satellite antenna model */
		//	if ( opt->posopt[0] )
		//	{
		satantpcv(rs + i * 6, rr, nav->pcvs + sat - 1, dants);
		//	}

		/* receiver antenna model */
		antmodel(opt->pcvr, opt->antdel[0], azel + i * 2, opt->posopt[1], dantr);

		/* phase windup correction */
		//	if ( opt->posopt[2] ) 
		//	{
		windupcorr(rtk->sol.time, rs + i * 6, rr, &rtk->ssat[sat - 1].phw);  //和频率无关
		//	}


		/* zenith hydrostatic delay */
		zhd = tropmodel(obs->time, pos, zazel, 0.0);  //humi==0 只计算干分量

		/* mapping function */
		m_h = tropmapf(obs->time, pos, azel + i * 2, &m_w);

		/* dcb correction  单频率如何进行dcb改正  */
		gamma = SQR(nav->lam[sat - 1][1] / nav->lam[sat - 1][0]); /* f1^2/f2^2 */

		P1_P2 = nav->cbias[sat - 1][0];
		P1_C1 = nav->cbias[sat - 1][1];

		clk_sat = CLIGHT*dts[i * 2];

		/* satellite clock and tropospheric delay ，这里需要扣除掉对流层的干分量*/
		r = geodis - clk_sat + (zhd*m_h);

		fprintf(g_mylog, "sat %02d ", sat);

		frequenCount = 0;
		//不用进行电离层改正，这里需要处理不同频率的硬件延迟的问题
		for (j = 0; j < opt->nf; j++)   //对于每个频率组建UofC观测值
		{
			test1 = (obs[i].L[j] - rtk->ssat[sat - 1].phw)*nav->lam[sat - 1][j] - dants[j] - dantr[j];  //载波
			test2 = obs[i].P[j] - dants[j] - dantr[j]; //伪距
			if (!(obs[i].L[j] != 0.0 && obs[i].P[j] != 0.0))
			{
				//f = myIB(sat,j,&rtk->opt);
				//rtk->x[f] = 0.0;  //模糊度设置为0
				//rtk->P[f+f*rtk->nx] =0.0;
				//f=0;
				continue;
			}

			v[nv] = (test1 + test2) / 2.0 - r;  // UofC, unit:m 观测值的方差跟伪距一样???,消电离层组合,硬件延迟比较复杂

			fprintf(g_mylog, "frequency %d ", j);
			fprintf(g_mylog, "%6.3f = L:%10.3f + C:%10.3f - G:%10.3f ,", v[nv], obs[i].L[j], obs[i].P[j], r);

			for (k = 0; k < nx; k++) H[k + nx*nv] = 0.0;  //先对矩阵H的nv行清0

			//方向余弦的问题(H)，实际上方向余弦存储在e数组中，但是还需要考虑动态情形
			for (k = 0; k < 3; k++) H[k + nx*nv] = -e[k];

			//然后扣除待估计参数的初值，滤波只估计参数的变化量

			v[nv] = v[nv] - x[IC(0, opt)];  //扣除接收机钟差的初始值			
			H[IC(0, opt) + nx*nv] = 1.0; //接收机钟差的设计矩阵系数;第nv行

			//扣除对流层参数(改正对流层干分量，只估计湿分量，因此这里应该是湿分量的投影系数 )
			k = IT(opt);
			v[nv] = v[nv] - x[IT(opt)] * m_w;  //这里的估计量是zwd
			H[IT(opt) + nx*nv] = m_w;   //湿分量投影系数

			//模糊度参数
			v[nv] = v[nv] - (x[myIB(obs[i].sat, j, opt)])*nav->lam[sat - 1][j] / 2.0;  //扣除掉模糊度的初始值；这里应该是模糊度的一半
			H[myIB(obs[i].sat, j, opt) + nx*nv] = (nav->lam[sat - 1][j]) / 2.0;  //模糊度的系数是波长的0.5倍

			var[nv] = 0.25*(varerr(sat, sys, azel[1 + i * 2], 1, opt) + varm[1] + vare[i] + vart);  //当前观测值的方差,这里的type强制设置为1,

			//设置卫星可用性
			rtk->ssat[sat - 1].vsat[j] = 1;

			nv++;

			frequenCount++;

		}  //end of for(j=0,...)

		if (frequenCount == 2)   //双频观测数据，需要再加一个载波之间的无电离层组合
		{
			if (!(obs[i].L[0] != 0.0 && obs[i].L[1] != 0.0))
			{
				continue;
			}

			alfa = gamma / (gamma - 1.0);
			beta = 1.0 / (gamma - 1.0);
			test1 = (obs[i].L[0] - rtk->ssat[sat - 1].phw)*nav->lam[sat - 1][0] - dants[0] - dantr[0];  // L1载波(m) 
			test2 = (obs[i].L[1] - rtk->ssat[sat - 1].phw)*nav->lam[sat - 1][1] - dants[1] - dantr[1];  // L2载波(m) 
			test1 = test1*alfa;
			test2 = test2*beta;
			v[nv] = test1 - test2 - r;   //载波无电离层组合观测值

			fprintf(g_mylog, "frequency %s ", "IF");
			fprintf(g_mylog, "%6.3f = %8.3f - %8.3f ,", v[nv], test1 - test2, r);

			for (k = 0; k < nx; k++) H[k + nx*nv] = 0.0;  //先对矩阵H的nv行清0

			//方向余弦的问题(H)，实际上方向余弦存储在e数组中，但是还需要考虑动态情形
			for (k = 0; k < 3; k++) H[k + nx*nv] = -e[k];

			//然后扣除待估计参数的初值，滤波只估计参数的变化量

			v[nv] = v[nv] - x[IC(0, opt)];  //扣除接收机钟差的初始值			
			H[IC(0, opt) + nx*nv] = 1.0; //无电离层组合接收机钟差系数仍然为1;第nv行

			//扣除对流层参数(改正对流层干分量，只估计湿分量，因此这里应该是湿分量的投影系数 )
			k = IT(opt);
			v[nv] = v[nv] - x[IT(opt)] * m_w;  //这里的估计量是zwd
			H[IT(opt) + nx*nv] = m_w;   //湿分量投影系数

			//模糊度参数
			v[nv] = v[nv] - (((x[myIB(obs[i].sat, 0, opt)])*nav->lam[sat - 1][0] * alfa)     //L1载波 
				- ((x[myIB(obs[i].sat, 1, opt)])*nav->lam[sat - 1][1] * beta)     //L2载波
				);  //扣除掉模糊度的初始值，注意这里全波长

			H[myIB(obs[i].sat, 0, opt) + nx*nv] = alfa*(nav->lam[sat - 1][0]);  //模糊度的系数是全波长*alfa
			H[myIB(obs[i].sat, 1, opt) + nx*nv] = -beta*(nav->lam[sat - 1][1]);  //模糊度的系数是全波长*beta

			//IF组合观测值方差
			var[nv] = (alfa*alfa + beta*beta)*
				(varerr(sat, sys, azel[1 + i * 2], 0, opt) + varm[1] + vare[i] + vart);  //当前观测值的方差,这里的type强制设置为0(载波),

			//设置卫星可用性
			rtk->ssat[sat - 1].vsat[j] = 1;

			nv++;
		}
		else if (opt->nf == 3)  //3频需要组合2个载波的无电离层组合
		{

		}

		fprintf(g_mylog, "\n");


	} // end of for(i=0,...)

	//设置观测值的方差协方差矩阵(各个频率不相关，因此是个方阵)
	for (i = 0; i < nv; i++)
	{
		for (j = 0; j < nv; j++)
		{
			R[i + j*nv] = i == j ? var[i] : 0.0;
		}
	}

	fprintf(g_mylog, "\n");

	return nv;
}




static obsd_t last_obs[MAXOBS * 2]; /* 保存上一历元的观测数据，用于进行周跳探测 */
static int    last_satnum;          //保存上一历元的卫星个数
static double *last_rs, *last_dts;   //保存上一历元的卫星位置及钟差
static double  tdobs[MAXSAT*NFREQ] = { 0.0 };  //每个频率都有的载波历元间差分观测值
static double  obsWeight[MAXSAT*NFREQ] = { 0.0 };
//rs = mat(6,satnum); dts=mat(2,satnum);

/*自定义的周跳探测函数
该周跳探测算法对于准静态效果很好
Author: Zhen.Li
*/
void mydetslp(int flag, rtk_t* rtk, const obsd_t* obs, const nav_t* nav, int satnum, double* rr, double* rs, double* dts)
{

	int maxIteratorCount = 10;  //最大迭代次数10次

	int i = 0, j = 0, f = 0, obscount = 0;
	int iterCount = 0;
	int maxsatnum = 0;
	double geodis_last = 0.0, geodis_now = 0.0, tropNow = 0.0, tropLast = 0.0;
	double e1[3] = { 0.0 }, e2[3] = { 0.0 }, pos[3] = { 0.0 };
	double rcvClkRate = 0.0, rcvClkRate_tmp = 0.0;
	double resdual[MAXSAT*NFREQ] = { 0.0 }, tmptotal = 0.0;  //残差
	double azel_now[2] = { 0.0 }, azel_last[2] = { 0.0 };
	double test1 = 0.0, test2 = 0.0, test = 0.0, test4 = 0.0, test5 = 0.0;

	maxsatnum = satnum >= last_satnum ? satnum : last_satnum;
	memset(tdobs, 0, sizeof(double)*MAXSAT*NFREQ);
	memset(obsWeight, 0, sizeof(double)*MAXSAT*NFREQ);

	ecef2pos(rr, pos);

	if (flag == 0)  //第一个历元不用进行周跳探测
	{
		last_satnum = satnum;
		memcpy(last_obs, obs, sizeof(obsd_t)*MAXOBS * 2);
		last_rs = mat(6, last_satnum);
		last_dts = mat(2, last_satnum);

		memcpy(last_rs, rs, sizeof(double)* 6 * satnum);
		memcpy(last_dts, dts, sizeof(double)* 2 * satnum);
	}
	else if (flag == 1)   //从第2个历元开始探测
	{

		for (i = 0; i < satnum; i++)  //当前历元所有卫星
		{
			geodis_now = geodist(rs + i * 6, rr, e1);  //当前历元测站坐标rr

			for (j = 0; j < last_satnum; j++)
			{
				geodis_last = geodist(last_rs + j * 6, rr, e2);  //上一历元测站坐标rr

				if (obs[i].sat == last_obs[j].sat)  //相同的卫星间作历元差分
				{
					satazel(pos, e1, azel_last);
					satazel(pos, e2, azel_now);

					//tropmodel(obs[i].time,)
					tropNow = tropmodel(obs[i].time, pos, azel_now, 0.7);
					tropLast = tropmodel(last_obs[j].time, pos, azel_last, 0.7);

					for (f = 0; f < rtk->opt.nf; f++)  //所有频率
					{
						if (fabs(obs[i].L[f]) < 100 || fabs(last_obs[j].L[f]) < 100)
						{
							continue;
						}

						test1 = (obs[i].L[f] - last_obs[j].L[f])*nav->lam[obs[i].sat - 1][f];  // uint: m
						test2 = (geodis_now - dts[i * 2] * CLIGHT  /*+ tropNow*/) - (geodis_last - last_dts[j * 2] * CLIGHT /* + tropLast*/);  //uint: m
						tdobs[(obs[i].sat - 1)*NFREQ + f] = test1 - test2;

						obscount++;
						obsWeight[(obs[i].sat - 1)*NFREQ + f] = 1.0;  //初始化时等权
						rcvClkRate = tdobs[(obs[i].sat - 1)*NFREQ + f] / obscount + rcvClkRate*(obscount - 1.0) / obscount; //等权均值

					}
					geodis_last = 0.0;
					geodis_now = 0.0;
					break;
				}
			}
		}





		//开始解算接收机钟差变化率(抗差估计)
		//本质上是加权平均值
		while (1)  //不断的调权处理
		{
			if (iterCount > maxIteratorCount)
			{
				break;
			}

			//计算残差并调权
			for (i = 0; i< MAXSAT*NFREQ; i++)
			{
				if (tdobs[i] != 0.0)
				{
					resdual[i] = tdobs[i] - rcvClkRate;
					if (fabs(resdual[i]) > 0.05)  //0残差大于.05m，则认为发生了周跳，开始调权
					{
						obsWeight[i] = 0.05 / resdual[i] / 100.0;
					}
				}
			}

			//对obsweight进行归一化处理
			for (i = 0; i < MAXSAT*NFREQ; i++)
			{
				tmptotal += obsWeight[i];
			}
			for (i = 0; i < MAXSAT*NFREQ; i++)
			{
				obsWeight[i] = obsWeight[i] / tmptotal;
			}

			//重新计算rcvClk
			for (i = 0; i < MAXSAT*NFREQ; i++)
			{
				rcvClkRate_tmp += (obsWeight[i] * tdobs[i]);
			}

			if (fabs(rcvClkRate_tmp - rcvClkRate) < 1.0E-4) //迭代结束
			{
				break;
			}
			else   //继续迭代
			{
				rcvClkRate = rcvClkRate_tmp;
			}

			iterCount++;
		}

		//对结果进行分析，找出每颗卫星每个频率的周跳值




		//最后更新数据
		last_satnum = satnum;
		memcpy(last_obs, obs, sizeof(obsd_t)*MAXOBS * 2);
		if (last_rs != NULL)  { free(last_rs); }
		if (last_dts != NULL)  { free(last_dts); }

		last_rs = mat(6, satnum);
		last_dts = mat(2, satnum);

		memcpy(last_rs, rs, sizeof(double)* 6 * satnum);
		memcpy(last_dts, dts, sizeof(double)* 2 * satnum);
	}





}

/*自定义的ppp数据处理核心函数
n实际上是卫星个数
Author: Zhen.Li
*/
extern int myPPP(rtk_t *rtk, const obsd_t *obs, int satnum, const nav_t *nav)
{
	prcopt_t *opt = &rtk->opt;
	sol_t solb = { { 0 } };
	gtime_t time = { 0 };
	int i, j, nu, nr, iTrop, iRcvClk;
	char msg[128] = "";
	char timestr[100] = "";

	double *rs, *dts, *var, *v, *H, *R, *azel, *xp, *Pp;
	int nv, info, svh[MAXOBS], stat = SOLQ_SINGLE;

	int is_first_epoch = 1;
	double test = 0.0;

	//	trace(3,"rtkpos  : time=%s n=%d\n",time_str(obs[0].time,3),n);
	//	trace(4,"obs=\n"); traceobs(4,obs,n);

	/* count rover/base station observations */
	//for ( nu=0;nu   <n&&obs[nu   ].rcv==1;nu++ );
	//for ( nr=0;nu+nr<n&&obs[nu+nr].rcv==2;nr++ );

	time = rtk->sol.time; /* previous epoch */

	time2str(obs[0].time, timestr, 100);
	fprintf(g_mylog, "time: %s\n", timestr);
	fprintf(g_mylog, "information: nx: %d \n", rtk->nx);

	if (norm(rtk->x, 3) <= 0.0)  //第一个历元进行
	{
		is_first_epoch = 0;

		//首先进行单点定位，这里必须采用广播星历
		/* rover position by single point positioning */
		if (!pntpos(obs, satnum, nav, &rtk->opt, &rtk->sol, NULL, rtk->ssat, msg))
		{
			printf("spp:pntpos error!\n");
		}
	}

	// time difference between  epochs
	if (time.time != 0) rtk->tt = timediff(rtk->sol.time, time);

	//下面直接进行PPP数据处理(单历元，多颗卫星)

	rs = mat(6, satnum); dts = mat(2, satnum); var = mat(1, satnum); azel = zeros(2, satnum);

	for (i = 0; i < MAXSAT; i++) rtk->ssat[i].fix[0] = 0;

	/* 估计参数排列顺序，前3个是坐标(XYZ)，接收机钟差，ZWD，以及模糊度*/

	//该函数需要改进
	//udstate_ppp(rtk,obs,n,nav);

	/* satellite positions and clocks with precise eph and clk*/
	satposs(obs[0].time, obs, satnum, nav, rtk->opt.sateph, rs, dts, var, svh);
	printf("epochTime: %s\n", timestr);
	for (i = 0; i < satnum; i++)
	{
		printf("sat:%02d X:%10.3f Y:%10.3f Z:%10.3f C:%10.3f\n", obs[i].sat, *(rs + 6 * i), *(rs + 6 * i + 1), *(rs + 6 * i + 2), *(dts + i * 2));
	}

	//这里记录下卫星轨道和钟差值
	//mydetslp( is_first_epoch, rtk, obs, nav, satnum, rtk->sol.rr, rs, dts );

	myUdstate_ppp(rtk, obs, satnum, nav);

	fprintf(g_mylog, "satellites(%2d): ", satnum);
	for (i = 0; i < satnum; i++)
	{
		fprintf(g_mylog, "%02d ", obs[i].sat);
	}
	fprintf(g_mylog, "\n");

	xp = mat(rtk->nx, 1);   //rtk->nx是所有待估计的参数个数 
	Pp = zeros(rtk->nx, rtk->nx);
	matcpy(xp, rtk->x, rtk->nx, 1); // rtk->x 是状态向量

	nv = satnum*(rtk->opt.nf) * 2;     // nv 是观测值个数(每个频率的载波和伪距观测值合成一个)  opt.nf 是频率个数
	v = mat(nv, 1);  // v 是观测值
	H = mat(rtk->nx, nv); // H是设计矩阵应该是nv行nx列
	R = mat(nv, nv);   //R是观测值方差协方差矩阵，应该是nv行nv列

	//首先进行非差的模型改正，然后再用改正过的观测值组成UofC观测值，然后再进行解算
	for (i = 0; i < 1; i++)  //迭代次数
	{
		nv = myPPP_RES(1, obs, satnum, rs, dts, var, svh, nav, xp, rtk, v, H, R, azel);
		/* measurement update */
		matcpy(Pp, rtk->P, rtk->nx, rtk->nx);

		fprintf(g_mylog, "initial for x:\n");
		for (j = 0; j < rtk->nx; j++)
		{
			if (rtk->x[j] != 0.0)
			{
				fprintf(g_mylog, "%10.3f ", rtk->x[j]);
			}
		}
		fprintf(g_mylog, "\n");

		if ((info = filter(xp, Pp, H, v, R, rtk->nx, nv)))
		{
			trace(2, "ppp filter error %s info=%d\n", time_str(rtk->sol.time, 0), info);
			break;
		}

		trace(4, "x(%d)=", i + 1); tracemat(4, xp, 1, NR(opt), 13, 4);

		stat = SOLQ_PPP;
	}

	if (stat == SOLQ_PPP)
	{
		/* postfit residuals */
		nv = myPPP_RES(1, obs, satnum, rs, dts, var, svh, nav, xp, rtk, v, H, R, azel);

		fprintf(g_mylog, "postfit residuals\n");
		for (i = 0; i < nv; i++)
		{
			fprintf(g_mylog, "%8.3f ", v[i]);
		}
		fprintf(g_mylog, "\n");

		/*
		update state and covariance matrix
		将滤波结果还原到x中
		*/
		matcpy(rtk->x, xp, rtk->nx, 1);
		matcpy(rtk->P, Pp, rtk->nx, rtk->nx);

		//这里进行ppp模糊度的固定
		printf("Solution: %.3f %.3f %.3f %.3f \n", rtk->x[0], rtk->x[1], rtk->x[2], rtk->x[3]);


		///* update solution status */
		rtk->sol.ns = 0;
		for (i = 0; i < satnum&&i < MAXOBS; i++)
		{
			/*if (!rtk->ssat[obs[i].sat-1].vsat[0])
			{
			continue;
			}*/
			rtk->ssat[obs[i].sat - 1].lock[0]++;
			rtk->ssat[obs[i].sat - 1].outc[0] = 0;
			rtk->ssat[obs[i].sat - 1].fix[0] = 4;
			rtk->sol.ns++;
		}

		rtk->sol.stat = stat;

		for (i = 0; i < 3; i++)  //坐标及其方差
		{
			rtk->sol.rr[i] = rtk->x[i];
			rtk->sol.qr[i] = (float)rtk->P[i + i*rtk->nx];
		}

		rtk->sol.dtr[0] = rtk->x[IC(0, opt)] / CLIGHT;  //接收机钟差(解算单位是m；存储单位是秒)

		rtk->sol.time = obs[0].time;

		//rtk->sol.ratio = (float)(rtk->P[IT(opt)*rtk->nx+IT(opt)]);
		rtk->sol.ratio = (float)(rtk->x[IT(opt)]);


		//方差
		rtk->sol.qr[0] = (float)rtk->P[0];  //c_xx
		rtk->sol.qr[1] = (float)rtk->P[rtk->nx + 1];  //c_yy

		rtk->sol.qr[2] = (float)rtk->P[2 * rtk->nx + 2];  //c_zz
		//协方差
		rtk->sol.qr[3] = (float)rtk->P[1];  //c_xy
		rtk->sol.qr[4] = (float)rtk->P[2 + rtk->nx]; //c_yz
		rtk->sol.qr[5] = (float)rtk->P[2];  //c_zx

		//rtk->sol.dtr[1]=rtk->x[IC(1,opt)]-rtk->x[IC(0,opt)];
		for (i = 0; i < satnum&&i < MAXOBS; i++)
		{
			rtk->ssat[obs[i].sat - 1].snr[0] = MIN(obs[i].SNR[0], obs[i].SNR[1]);
		}
		for (i = 0; i < MAXSAT; i++)
		{
			if (rtk->ssat[i].slip[0] & 3)
			{
				rtk->ssat[i].slipc[0]++;
			}
		}
	}

	return 0;

}