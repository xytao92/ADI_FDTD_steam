//writen by liuzhichao 
//2016/10/8
#ifndef _GRID_
#define _GRID_

struct Grid	
{
    double ex,ey,ez;//电场强度
	double bx,by,bz;//磁场强度	
};
//对象统一设为一维的数组，方便后边并行化

Grid* halfgrid_now = new  Grid[Nx*Ny*Nz]; //半整数网格点现在时刻
Grid* halfgrid_before = new Grid[Nx*Ny*Nz];//半整数网格点过去时刻
Grid* halfgrid_beforeX2 = new  Grid[Nx*Ny*Nz];//半整数网格点过去的过去时刻,mur吸收边界要用

//Grid* grid_result = new Grid[Nx*Ny*Nz];//最终网格中场的结果

//Grid halfgrid_now[Ny*Nz*Nz]; //半整数网格点现在时刻
//Grid halfgrid_before[Ny*Nz*Nz];//半整数网格点过去时刻
//Grid grid_result[Ny*Nz*Nz];//最终网格中场结果

#endif