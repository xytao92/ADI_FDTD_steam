/*
Writen By LiuZhiChao
on time: 2017.3.14
All Rights Receved
*/

#ifndef _GSS_ADI_CPPML_
#define _GSS_ADI_CPPML_

#include "definer.h"
#include "grid.h"
#include "functions.h"
#include "Gss-2.0.h"
//------------------------------------------函数声明
void adi_fdtd_leapforg_cpml_GSS(Grid* halfgrid_now);
void gss_cal_cpml(Grid* halfgrid_now, int step);

void norm_gsscalc_ez(Grid* halfgrid_now, int step);
void norm_gsscalc_ex(Grid* halfgrid_now, int step);
void norm_gsscalc_ey(Grid* halfgrid_now, int step);

void norm_gsscalc_bz(Grid* halfgrid_now, int step);
void norm_gsscalc_bx(Grid* halfgrid_now, int step);
void norm_gsscalc_by(Grid* halfgrid_now, int step);

void cpml_gsscalc_ez(Grid* halfgrid_now, int step);
void cpml_gsscalc_ex(Grid* halfgrid_now, int step);
void cpml_gsscalc_ey(Grid* halfgrid_now, int step);

void cpml_gsscalc_bz(Grid* halfgrid_now, int step);
void cpml_gsscalc_bx(Grid* halfgrid_now, int step);
void cpml_gsscalc_by(Grid* halfgrid_now, int step);


//-------------常量声明
//CPML常数部分
const int kxmax = 8;
const int kymax = 8;
const int kzmax = 8;

const double sigma_xmax = (4 + 1) / (150 * pi*dz);
const double alpha_xmax = 0.05;

const double sigma_ymax = (4 + 1) / (150 * pi*dz);
const double alpha_ymax = 0.05;

const double sigma_zmax = (4 + 1) / (150 * pi*dz);
const double alpha_zmax = 0.05;

const double d = 10 * dz;//cpml层的总厚度

#endif