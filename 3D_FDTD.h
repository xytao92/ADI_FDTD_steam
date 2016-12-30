#ifndef _3D_FDTD_
#define _3D_FDTD_

#include "definer.h"
#include "grid.h"
#include"functions.h"
#include"Gss-2.0.h"

//-----------------------------------------------------------------------函数声明
void fdtd_ez(Grid* halfgrid_before, Grid* halfgrid_now, int step);
void fdtd_ex(Grid* halfgrid_before, Grid* halfgrid_now, int step);
void fdtd_ey(Grid* halfgrid_before, Grid* halfgrid_now, int step);

void fdtd_bz(Grid* halfgrid_before, Grid* halfgrid_now, int step);
void fdtd_bx(Grid* halfgrid_before, Grid* halfgrid_now, int step);
void fdtd_by(Grid* halfgrid_before, Grid* halfgrid_now, int step);

void fdtd_matel(Grid* halfgrid_before, Grid* halfgrid_now)

{
	//system("mkdir result");

	ofstream file("result\\fdtd_01.txt");//
	ofstream file2("result\\fdtd_plat_01.txt");//

	//*******计算TE10模******//注意边界条件的问题，不处理周围的四个面												   
	int result_z = 2;
	int result_x = 10;
	int result_y = 5;

	int step = 0;//计算时间步长
	inject_field(halfgrid_before, halfgrid_now, step);

	while (step < STEPS)
	{		
//---------------------------------------------------------------------------------------计算磁场
		fdtd_bz(halfgrid_before, halfgrid_now, step);
		fdtd_bx(halfgrid_before, halfgrid_now, step);
		fdtd_by(halfgrid_before, halfgrid_now, step);
//---------------------------------------------------------------------------------------计算电场
		fdtd_ez(halfgrid_before, halfgrid_now, step);
		fdtd_ex(halfgrid_before, halfgrid_now, step);
		fdtd_ey(halfgrid_before, halfgrid_now, step);

		file << step << '\t' << halfgrid_now[result_z * Nx*Ny + result_x * Ny + result_y].ex 
			 << '\t' << halfgrid_now[result_z * Nx*Ny + result_x * Ny + result_y].ey << '\t' 
			 << halfgrid_now[result_z * Nx*Ny + result_x * Ny + result_y].ez << '\t';

		file << halfgrid_now[result_z * Nx*Ny + result_x * Ny + result_y].bx 
			<< '\t' << halfgrid_now[result_z * Nx*Ny + result_x * Ny + result_y].by 
			<< '\t' << halfgrid_now[result_z * Nx*Ny + result_x * Ny + result_y].bz << '\t';

		file << '\n';

		cout << "Step--- " << step << " ---has finished." << endl;
		step++;
	}//while

	for (int i = 0; i < Nx - 1; i++)//输出一个横截面的数据
	{
		for (int k = 0; k < Nz - 1; k++)
		{
			file2 << k << '\t' << i << '\t' 
				  << halfgrid_now[k*Nx*Ny + i*Ny + 5].ey << '\t' 
				  << halfgrid_now[k*Nx*Ny + i*Ny + 5].bx << '\t' 
				  << halfgrid_now[k*Nx*Ny + i*Ny + 5].bz << endl;
		}

	}

}//函数结尾

void fdtd_ez(Grid* halfgrid_before, Grid* halfgrid_now, int step)
{
	for (int i = 0; i < Nx - 2; i++)
	{
		for (int j = 0; j < Ny - 2; j++)
		{
			for (int k = 0; k < Nz - 2; k++)
			{
				halfgrid_before[k * Nx*Ny + i * Ny + j].ez = halfgrid_now[k * Nx*Ny + i * Ny + j].ez;
				if (j == 0 || j == Ny - 2 || i == 0 || i == Nx - 2)
				{
					halfgrid_now[k * Nx*Ny + i * Ny + j].ez = 0.0;
				}
				else
				{
					halfgrid_now[k * Nx*Ny + i * Ny + j].ez = halfgrid_before[k * Nx*Ny + i * Ny + j].ez
						+ (dt / epsl0)*
						((halfgrid_now[k * Nx*Ny + i * Ny + j].by - halfgrid_now[k * Nx*Ny + (i - 1) * Ny + j].by) / dx
							- (halfgrid_now[k * Nx*Ny + i * Ny + j].bx - halfgrid_now[k * Nx*Ny + i * Ny + j - 1].bx) / dy);
				}
				
			}//for
		}//for
	}//for
	
}

void fdtd_ex(Grid* halfgrid_before, Grid* halfgrid_now, int step)
{
	for (int i = 0; i < Nx - 2; i++)
	{
		for (int j = 0; j < Ny - 2; j++)
		{
			for (int k = 0; k < Nz - 2; k++)
			{
				halfgrid_before[k * Nx*Ny + i * Ny + j].ex = halfgrid_now[k * Nx*Ny + i * Ny + j].ex;
				if (j == 0 || j == Ny - 2 || k == Nz - 2 || k == 0)
				{
					halfgrid_now[k * Nx*Ny + i * Ny + j].ex = 0.0;
				}
				else
				{
					halfgrid_now[k * Nx*Ny + i * Ny + j].ex = halfgrid_before[k * Nx*Ny + i * Ny + j].ex
						+ (dt / epsl0)*
						((halfgrid_now[k * Nx*Ny + i * Ny + j].bz - halfgrid_now[k * Nx*Ny + i * Ny + j-1].bz) / dy
							- (halfgrid_now[k * Nx*Ny + i * Ny + j].by - halfgrid_now[(k-1) * Nx*Ny + i * Ny + j].by) / dz);
				}

			}//for
		}//for
	}//for

		
}


void fdtd_ey(Grid* halfgrid_before, Grid* halfgrid_now, int step)
{
	for (int i = 0; i < Nx - 2; i++)
	{
		for (int j = 0; j < Ny - 2; j++)
		{
			for (int k = 0; k < Nz - 2; k++)
			{
				halfgrid_before[k * Nx*Ny + i * Ny + j].ey = halfgrid_now[k * Nx*Ny + i * Ny + j].ey;
				if (k == 0)
				{
					if (step*dt < 2 * T)
					{
						halfgrid_now[k * Nx*Ny + i * Ny + j].ey = ((step*dt) / (2 * T))*((omega*mur0*X) / pi) * hm * sin((pi / X)*i*dx) * sin(omega*step*dt);
					}
					else
					{
						halfgrid_now[k * Nx*Ny + i * Ny + j].ey = ((omega*mur0*X) / pi) * hm * sin((pi / X)*i*dx)*sin(omega*step*dt);//不改变源的值
					}
				}

				if (i == 0 || i == Nx - 2 || k == Nz - 2)
				{
					halfgrid_now[k * Nx*Ny + i * Ny + j].ey = 0.0;
				}
				else
				{
					halfgrid_now[k * Nx*Ny + i * Ny + j].ey = halfgrid_before[k * Nx*Ny + i * Ny + j].ey
						+ (dt / epsl0)*
						((halfgrid_now[k * Nx*Ny + i * Ny + j].bx - halfgrid_now[(k-1) * Nx*Ny + i * Ny + j].bx) / dz
							- (halfgrid_now[k * Nx*Ny + i * Ny + j].bz - halfgrid_now[k * Nx*Ny + (i-1) * Ny + j].bz) / dx);
				}

			}//for
		}//for
	}//for
			
}


void fdtd_bz(Grid* halfgrid_before, Grid* halfgrid_now, int step)
{
	for (int i = 0; i < Nx - 2; i++)
	{
		for (int j = 0; j < Ny - 2; j++)
		{
			for (int k = 0; k < Nz - 2; k++)
			{
				halfgrid_before[k * Nx*Ny + i * Ny + j].bz = halfgrid_now[k * Nx*Ny + i * Ny + j].bz;
				if (k == 0)//不改变输入源的值
				{
					if (step*dt < 2 * T)
					{
						halfgrid_now[k * Nx*Ny + i * Ny + j].ez = ((step*dt) / (2 * T))*hm*cos((pi / X)*i*dx)*cos(omega*step*dt*1.5);
					}
					else
					{
						halfgrid_now[k * Nx*Ny + i * Ny + j].ez = hm * cos((pi / X)*i*dx)*cos(omega*step*dt*1.5);
					}

				}
				if (k==Nz-2)
				{
					halfgrid_now[k * Nx*Ny + i * Ny + j].bz = 0.0;
				}
				else
				{
					halfgrid_now[k * Nx*Ny + i * Ny + j].bz = halfgrid_before[k * Nx*Ny + i * Ny + j].bz
						+ (dt / mur0)*
						((halfgrid_before[k * Nx*Ny + (i+1) * Ny + j].ey - halfgrid_before[k * Nx*Ny + i * Ny + j].ey) / dx
							- (halfgrid_before[k * Nx*Ny + i * Ny + j+1].ex - halfgrid_before[k * Nx*Ny + i * Ny + j].ex) / dy);
				}

			}//for
		}//for
	}//for
			
}

void fdtd_bx(Grid* halfgrid_before, Grid* halfgrid_now, int step)
{
	for (int i = 0; i < Nx - 2; i++)
	{
		for (int j = 0; j < Ny - 2; j++)
		{
			for (int k = 0; k < Nz - 2; k++)
			{	
				halfgrid_before[k * Nx*Ny + i * Ny + j].bx = halfgrid_now[k * Nx*Ny + i * Ny + j].bx;

					halfgrid_now[k * Nx*Ny + i * Ny + j].bx = halfgrid_before[k * Nx*Ny + i * Ny + j].bx
						+ (dt / epsl0)*
						((halfgrid_before[k * Nx*Ny + i * Ny + j+1].ez - halfgrid_before[k * Nx*Ny + i * Ny + j].ez) / dy
							- (halfgrid_before[(k+1) * Nx*Ny + i * Ny + j].ey - halfgrid_before[k * Nx*Ny + i * Ny + j].ey) / dz);
				
			}//for
		}//for
	}//for
			
}

void fdtd_by(Grid* halfgrid_before, Grid* halfgrid_now, int step)
{
	for (int i = 0; i < Nx - 2; i++)
	{
		for (int j = 0; j < Ny - 2; j++)
		{
			for (int k = 0; k < Nz - 2; k++)
			{
				halfgrid_before[k * Nx*Ny + i * Ny + j].by = halfgrid_now[k * Nx*Ny + i * Ny + j].by;

					halfgrid_now[k * Nx*Ny + i * Ny + j].by = halfgrid_before[k * Nx*Ny + i * Ny + j].by
						+ (dt / epsl0)*
						((halfgrid_before[(k+1) * Nx*Ny + i * Ny + j].ex - halfgrid_before[k * Nx*Ny + i * Ny + j].ex) / dz
							- (halfgrid_before[k * Nx*Ny + (i+1) * Ny + j].ez - halfgrid_before[k * Nx*Ny + i * Ny + j].ez) / dx);
				
			}//for
		}//for
	}//for

}




#endif