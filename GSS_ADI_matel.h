
/*
Writen By LiuZhiChao
on time: 2016.11.22
All Rights Receved
*/

#ifndef _GSS_ADI_MATEL_
#define _GSS_ADI_MATEL_

#include "definer.h"
#include "grid.h"
#include "functions.h"
#include "Gss-2.0.h"
//------------------------------------------º¯ÊýÉùÃ÷
void adi_fdtd_leapforg_matel_GSS(Grid* halfgrid_now);
void gss_cal_matel(Grid* halfgrid_now, int step);

void matel_gsscalc_ez(Grid* halfgrid_now, int step);
void matel_gsscalc_ex(Grid* halfgrid_now, int step);
void matel_gsscalc_ey(Grid* halfgrid_now, int step);

void matel_gsscalc_bz(Grid* halfgrid_now, int step);
void matel_gsscalc_bx(Grid* halfgrid_now, int step);
void matel_gsscalc_by(Grid* halfgrid_now, int step);








#endif
