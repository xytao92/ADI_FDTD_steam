#ifndef _3D_FDTD_
#define _3D_FDTD_

#include "definer.h"
#include "grid.h"
#include"functions.h"
#include"Gss-2.0.h"
//1.3上午检查了公式的问题并做了修正，修正了循环范围的问题
//1.4上午继续调试FDTD代码，检查满足CFL稳定条件？->满足。存储过程？->发现根本不需要halfgrid-before.直接用now更新即可
//边界处理？循环控制？参数设置？
//-----------------------------------------------------------------------函数声明
void fdtd_bz(Grid* halfgrid_now, int step);
void fdtd_bx(Grid* halfgrid_now, int step);
void fdtd_by(Grid* halfgrid_now, int step);

void fdtd_ez(Grid* halfgrid_now, int step);
void fdtd_ex(Grid* halfgrid_now, int step);
void fdtd_ey(Grid* halfgrid_now, int step);

void fdtd_matel(Grid* halfgrid_now);
void fdtd_bz(Grid* halfgrid_now, int step);


#endif