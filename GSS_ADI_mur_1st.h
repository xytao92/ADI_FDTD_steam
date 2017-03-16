
/*
Writen By LiuZhiChao
on time: 2016.11.22
All Rights Receved
*/

#ifndef _GSS_ADI_MUR_1st_
#define _GSS_ADI_MUR_1st_

#include "definer.h"
#include "grid.h"
#include"functions.h"
#include"Gss-2.0.h"

//------------------------------------------函数声明
void adi_fdtd_leapforg_mur1_GSS(Grid* halfgrid_now);
void gss_cal_mur1(Grid* halfgrid_now, int step);

void mur1_gsscalc_ez(Grid* halfgrid_now, int step);
void mur1_gsscalc_ex(Grid* halfgrid_now, int step);
void mur1_gsscalc_ey(Grid* halfgrid_now, int step);

void mur1_gsscalc_bz(Grid* halfgrid_now, int step);
void mur1_gsscalc_bx(Grid* halfgrid_now, int step);
void mur1_gsscalc_by(Grid* halfgrid_now, int step);





//=======================输出横截面数据========================//

#endif
