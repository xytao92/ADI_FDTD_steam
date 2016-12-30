//writen by liuzhichao 
//2016/10/8

#ifndef _FUNCTIONS_
#define _FUNCTIONS_
#include "definer.h"
#include "grid.h"


void initGrid(Grid* halfgrid_beforeX2,Grid* halfgrid_before, Grid* halfgrid_now)//网格的初始化
{
	for (int k = 0; k<Nz; k++)
		for (int i = 0; i<Nx; i++)
			for (int j = 0; j<Ny; j++)
			{
				halfgrid_now[k*Nx*Ny + i*Ny + j].ex = 0.0;
				halfgrid_now[k*Nx*Ny + i*Ny + j].ey = 0.0;
				halfgrid_now[k*Nx*Ny + i*Ny + j].ez = 0.0;
				halfgrid_now[k*Nx*Ny + i*Ny + j].bx = 0.0;
				halfgrid_now[k*Nx*Ny + i*Ny + j].by = 0.0;
				halfgrid_now[k*Nx*Ny + i*Ny + j].bz = 0.0;

				halfgrid_before[k*Nx*Ny + i*Ny + j].ex = 0.0;
				halfgrid_before[k*Nx*Ny + i*Ny + j].ey = 0.0;
				halfgrid_before[k*Nx*Ny + i*Ny + j].ez = 0.0;
				halfgrid_before[k*Nx*Ny + i*Ny + j].bx = 0.0;
				halfgrid_before[k*Nx*Ny + i*Ny + j].by = 0.0;
				halfgrid_before[k*Nx*Ny + i*Ny + j].bz = 0.0;

				halfgrid_beforeX2[k*Nx*Ny + i*Ny + j].ex = 0.0;
				halfgrid_beforeX2[k*Nx*Ny + i*Ny + j].ey = 0.0;
				halfgrid_beforeX2[k*Nx*Ny + i*Ny + j].ez = 0.0;
				halfgrid_beforeX2[k*Nx*Ny + i*Ny + j].bx = 0.0;
				halfgrid_beforeX2[k*Nx*Ny + i*Ny + j].by = 0.0;
				halfgrid_beforeX2[k*Nx*Ny + i*Ny + j].bz = 0;

				/*grid_result[k*Nx*Ny + i*Ny + j].ex = 0.0;
				grid_result[k*Nx*Ny + i*Ny + j].ey = 0.0;
				grid_result[k*Nx*Ny + i*Ny + j].ez = 0.0;
				grid_result[k*Nx*Ny + i*Ny + j].bx = 0.0;
				grid_result[k*Nx*Ny + i*Ny + j].by = 0.0;
				grid_result[k*Nx*Ny + i*Ny + j].bz = 0.0;*/
			}

}

void free(Grid* halfgrid_beforeX2,Grid* halfgrid_before, Grid* halfgrid_now)
{
	delete halfgrid_before;
	delete halfgrid_now;
	delete halfgrid_beforeX2;
	//delete grid_result;
}

void inject_field(Grid* halfgrid_before, Grid* halfgrid_now,int step)//计算激励源//
{
	int k =0;//加点源的位置//
	
	//halfgrid_before[k*Ny*Nz + i*Nz + j].ez = hm*sin(omega*step*dt);//hm为设置值//115
	//halfgrid_now[k*Nx*Ny + i*Ny + j].ez = hm * sin(omega*step*dt);
	
	for (int i = 1; i < Nx - 2; i++)
	{
		for (int j = 1; j < Ny - 2; j++)//不能将源加到边界面上
		{
			if (step*dt < 2 * T)
			{
				//加一个完整的TE10,模式 电场-----------------------------------------------------
				/*halfgrid_beforeX2[i*Ny + j].ez = ((step*dt) / (2 * T))*0.0;
				halfgrid_before[i*Ny + j].ez =halfgrid_beforeX2[i*Ny + j].ez;
				halfgrid_now[i*Ny + j].ez = halfgrid_before[i*Ny + j].ez;*/

				halfgrid_beforeX2[i*Ny + j].ex = ((step*dt) / (2 * T))*0.0;
				halfgrid_before[i*Ny + j].ex = halfgrid_beforeX2[i*Ny + j].ex;
				halfgrid_now[i*Ny + j].ex = halfgrid_before[i*Ny + j].ex;

				halfgrid_beforeX2[i*Ny + j].ey = ((step*dt) / (2 * T))*((omega*mur0*X) / pi) * hm * sin((pi / X)*i*dx) * sin(omega*step*dt); //-1 * pi*sin(pi*i*dx / X)*sin(2 * pi*freq*step*dt) / X;//((omega*mur0*X) / pi) * hm * sin((pi / X)*i*dx) * sin(omega*step*dt);
				halfgrid_before[i*Ny + j].ey = halfgrid_beforeX2[i*Ny + j].ey;//hm为设置值//115
				halfgrid_now[i*Ny + j].ey = halfgrid_before[i*Ny + j].ey;//hm为设置值//115

																		 //磁场---------------------------------------------------------------------------
				halfgrid_beforeX2[i*Ny + j].bz = ((step*dt) / (2 * T))*hm*cos((pi / X)*i*dx)*cos(omega*step*dt*1.5);
				halfgrid_before[i*Ny + j].bz = halfgrid_beforeX2[i*Ny + j].bz;
				halfgrid_now[i*Ny + j].bz = halfgrid_before[i*Ny + j].bz;


				/*halfgrid_beforeX2[i*Ny + j].bx = -1 *((step*dt) / (2 * T))* (X*bate / pi)*hm*sin((pi / X)*i*dx)*sin(omega*step*dt*1.5);
				halfgrid_before[i*Ny + j].bx = halfgrid_beforeX2[i*Ny + j].bx;
				halfgrid_now[i*Ny + j].bx = halfgrid_before[i*Ny + j].bx;*/

				/*halfgrid_beforeX2[i*Ny + j].bx =((step*dt) / (2 * T))*0.0;
				halfgrid_before[i*Ny + j].bx = halfgrid_beforeX2[i*Ny + j].bx;
				halfgrid_now[i*Ny + j].bx = halfgrid_before[i*Ny + j].bx;*/

				/*halfgrid_beforeX2[i*Ny + j].by = ((step*dt) / (2 * T))*0.0;
				halfgrid_before[i*Ny + j].by = halfgrid_beforeX2[i*Ny + j].by;
				halfgrid_now[i*Ny + j].by = halfgrid_before[i*Ny + j].by;*/

			}
			else 
			{
				//加一个完整的TE10,模式 电场-----------------------------------------------------
				/*halfgrid_beforeX2[i*Ny + j].ez = 0.0;
				halfgrid_before[i*Ny + j].ez =halfgrid_beforeX2[i*Ny + j].ez;
				halfgrid_now[i*Ny + j].ez = halfgrid_before[i*Ny + j].ez;*/

				halfgrid_beforeX2[i*Ny + j].ex = 0.0;
				halfgrid_before[i*Ny + j].ex = halfgrid_beforeX2[i*Ny + j].ex;
				halfgrid_now[i*Ny + j].ex = halfgrid_before[i*Ny + j].ex;

				halfgrid_beforeX2[i*Ny + j].ey = ((omega*mur0*X) / pi) * hm * sin((pi / X)*i*dx) * sin(omega*step*dt); //-1 * pi*sin(pi*i*dx / X)*sin(2 * pi*freq*step*dt) / X;//((omega*mur0*X) / pi) * hm * sin((pi / X)*i*dx) * sin(omega*step*dt);
				halfgrid_before[i*Ny + j].ey = halfgrid_beforeX2[i*Ny + j].ey;//hm为设置值//115
				halfgrid_now[i*Ny + j].ey = halfgrid_before[i*Ny + j].ey;//hm为设置值//115

																		 //磁场---------------------------------------------------------------------------
				halfgrid_beforeX2[i*Ny + j].bz = hm*cos((pi / X)*i*dx)*cos(omega*step*dt*1.5);
				halfgrid_before[i*Ny + j].bz = halfgrid_beforeX2[i*Ny + j].bz;
				halfgrid_now[i*Ny + j].bz = halfgrid_before[i*Ny + j].bz;


				/*halfgrid_beforeX2[i*Ny + j].bx = -1 * (X*bate / pi)*hm*sin((pi / X)*i*dx)*sin(omega*step*dt*1.5);
				halfgrid_before[i*Ny + j].bx = halfgrid_beforeX2[i*Ny + j].bx;
				halfgrid_now[i*Ny + j].bx = halfgrid_before[i*Ny + j].bx;*/

				/*halfgrid_beforeX2[i*Ny + j].bx =0.0;
				halfgrid_before[i*Ny + j].bx = halfgrid_beforeX2[i*Ny + j].bx;
				halfgrid_now[i*Ny + j].bx = halfgrid_before[i*Ny + j].bx;*/

				/*halfgrid_beforeX2[i*Ny + j].by = 0.0;
				halfgrid_before[i*Ny + j].by = halfgrid_beforeX2[i*Ny + j].by;
				halfgrid_now[i*Ny + j].by = halfgrid_before[i*Ny + j].by;*/
			
			}
			
		}
	}
}



//以下函数没有修改坐标系，按照老的坐标系来------------------------！！！！！！！

void adi_fdtd_leapforg_matel(Grid* halfgrid_before,Grid* halfgrid_now)

{
	/*system("mkdir result");*/
	ofstream file("result\\leapforg_ADI_FDTD_steam0.1_matel03.txt");//用于保存结果
	int i1 = 0;
	int j1 = 0;
	int k1 = 0;
	//*******计算TE10模******//注意边界条件的问题，不处理周围的四个面
	//PART1---- 计算电场//

	int step = 0;//计算时间步长
	while(step< STEPS)
	{
		if (step == 50)
		{
			cout << "It's Time, step =50!" << endl;
		}

		inject_field(halfgrid_before,halfgrid_now,step);
		for(int i = 0; i<Nx-1;i++)
			for(int j =1; j<Ny-2; j++)
				for (int k = 1; k<Nz-2; k++)
				{
			i1 = i;
			j1 = j;
			k1 = k;
			if (i == 100 && j == 100 && k == 50)
			{
				cout << "Now" << endl;
			}
			if (i == 0)//左边界，信号源处
			{
				i1 = 1;
			}
			
			//ex
			double ty0 = ( halfgrid_before[i*Ny*Nz + j*Nz + k].bz - halfgrid_before[i*Ny*Nz + (j1-1)*Nz + k].bz )/dy;
			double tz0 = ( halfgrid_before[i*Ny*Nz + j*Nz + k].by-halfgrid_before[i*Ny*Nz + j*Nz + k1-1].by )/dz;
			halfgrid_now[ i*Ny*Nz + j*Nz + k].ex = halfgrid_before[i*Ny*Nz + j*Nz + k].ex + dt*(1/epsl_x)*( ty0 - tz0 );
			//ey
			double tx1 =  ( halfgrid_before[i*Ny*Nz + j*Nz + k].bx - halfgrid_before[i*Ny*Nz + j*Nz + k1-1].bx )/dz;
			double tz1 =  ( halfgrid_before[i*Ny*Nz + j*Nz + k].bz -halfgrid_before[(i1-1)*Ny*Nz + j*Nz + k].bz )/dx;
			halfgrid_now[ i*Ny*Nz + j*Nz + k].ey = halfgrid_before[ i*Ny*Nz + j*Nz + k].ey + dt*(1/epsl_y)*( tx1 - tz1 );
			//ez
			double ty2 =  ( halfgrid_before[i*Ny*Nz + j*Nz + k].by - halfgrid_before[(i1-1)*Ny*Nz + j*Nz + k ].by )/dx;
			double tx2 =  ( halfgrid_before[i*Ny*Nz + j*Nz + k].bx-halfgrid_before[i*Ny*Nz + (j1-1)*Nz + k].bx )/dy;
			halfgrid_now[ i*Ny*Nz + j*Nz + k].ez = halfgrid_before[ i*Ny*Nz + j*Nz + k].ez + dt*(1/epsl_z)*( ty2 - tx2 );

			halfgrid_before[ i*Ny*Nz + j*Nz + k].ex = halfgrid_now[ i*Ny*Nz + j*Nz + k].ex;//保存为前一步结果
			halfgrid_before[ i*Ny*Nz + j*Nz + k].ey =  halfgrid_now[i*Ny*Nz + j*Nz + k].ey;
			halfgrid_before[i*Ny*Nz + j*Nz + k].ez  = halfgrid_now[i*Ny*Nz + j*Nz + k].ez;

			//bx
			double ty3 = (halfgrid_before[i*Ny*Nz + j*Nz + k + 1].ey - halfgrid_before[i*Ny*Nz + j*Nz + k].ey) / dz;
			double tz3 = (halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k].ez - halfgrid_before[i*Ny*Nz + j*Nz + k].ez) / dy;
			halfgrid_now[i*Ny*Nz + j*Nz + k].bx = halfgrid_before[i*Ny*Nz + j*Nz + k].bx + dt*(1 / mur_x)*(ty3 - tz3);
			//by
			double tz4 = (halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k].ez - halfgrid_before[i*Ny*Nz + j*Nz + k].ez) / dx;
			double tx4 = (halfgrid_before[i*Ny*Nz + j*Nz + k + 1].ex - halfgrid_before[i*Ny*Nz + j*Nz + k].ex) / dz;
			halfgrid_now[i*Ny*Nz + j*Nz + k].by = halfgrid_before[i*Ny*Nz + j*Nz + k].by + dt*(1 / mur_y)*(tz4 - tx4);
			//bz
			double tx5 = (halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k].ex - halfgrid_before[i*Ny*Nz + j*Nz + k].ex) / dy;
			double ty5 = (halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k].ey - halfgrid_before[i*Ny*Nz + j*Nz + k].ey) / dx;
			halfgrid_now[i*Ny*Nz + j*Nz + k].bz = halfgrid_before[i*Ny*Nz + j*Nz + k].bz + dt*(1 / mur_z)*(tx5 - ty5);

			halfgrid_before[i*Ny*Nz + j*Nz + k].bx = halfgrid_now[i*Ny*Nz + j*Nz + k].bx;//保存为前一步结果
			halfgrid_before[i*Ny*Nz + j*Nz + k].by = halfgrid_now[i*Ny*Nz + j*Nz + k].by;
			halfgrid_before[i*Ny*Nz + j*Nz + k].bz = halfgrid_now[i*Ny*Nz + j*Nz + k].bz;
			
			
//PART2----计算磁场//
	//for(int i = 0; i<Nx-1;i++)
	//	for(int j =0; j<Ny-1; j++)
	//		for (int k = 0; k < Nz-1; k++)
	//		{
	//			//bx
	//			double ty3 =  ( halfgrid_before[i*Ny*Nz + j*Nz + k+1].ey - halfgrid_before[ i*Ny*Nz + j*Nz + k].ey )/dz;
	//			double tz3 =  ( halfgrid_before[ i*Ny*Nz + (j+1)*Nz + k].ez-halfgrid_before[ i*Ny*Nz +  j*Nz + k ].ez  )/dy;
	//			halfgrid_now[i*Ny*Nz + j*Nz + k].bx = halfgrid_before[ i*Ny*Nz + j*Nz + k].bx + dt*(1/ mur_x )*( ty3 - tz3 );
	//			//by
	//			double tz4 =  ( halfgrid_before[ (i+1)*Ny*Nz + j*Nz + k].ez - halfgrid_before[ i*Ny*Nz +  j*Nz + k].ez )/dx;
	//			double tx4 =  (halfgrid_before[i*Ny*Nz + j*Nz + k+1].ex - halfgrid_before[i*Ny*Nz + j*Nz + k].ex )/dz;
	//			halfgrid_now[i*Ny*Nz + j*Nz + k].by = halfgrid_before[i*Ny*Nz + j*Nz + k].by + dt*(1/ mur_y )*( tz4 - tx4 );
	//			//bz
	//			double tx5 =  ( halfgrid_before[ i*Ny*Nz + (j+1)*Nz + k].ex - halfgrid_before[ i*Ny*Nz + j*Nz + k ].ex )/dy;
	//			double ty5 =  (halfgrid_before[ (i+1)*Ny*Nz + j*Nz + k].ey - halfgrid_before[ i*Ny*Nz + j*Nz + k].ey )/dx;
	//			halfgrid_now[ i*Ny*Nz + j*Nz + k].bz = halfgrid_before[ i*Ny*Nz + j*Nz + k].bz + dt*(1/ mur_z )*( tx5 - ty5 );			

	//			halfgrid_before[i*Ny*Nz + j*Nz + k].bx = halfgrid_now[i*Ny*Nz + j*Nz + k].bx;//保存为前一步结果
 //               halfgrid_before[ i*Ny*Nz + j*Nz + k].by =  halfgrid_now[i*Ny*Nz + j*Nz + k].by;
	//			halfgrid_before[ i*Ny*Nz + j*Nz + k].bz  = halfgrid_now[ i*Ny*Nz + j*Nz + k].bz;
			}
								
	////若取平均，则网格点处的电磁场计算如下
	//
	//for(int i = 0; i<Nx-1;i++)
	//	for(int j =0; j<Ny-1; j++)
	//		for (int k = 0; k < Nz-1; k++)
	//		{
	//			  if(i == 0 )//左边界，信号源处
	//			 {
	//				
	//			 }
	//		           
	//		     if(i== Nx - 1 )//右边界，吸收边界
	//			 {
	//				
	//			 }
	//			 if(k == Nz - 1|| k==0 || j==Ny - 1 || j==0 )//四个平面
	//			 {
	//				 
	//			 }
	//			 else
	//			 {
	//			 
	//			    grid_result[i*Ny*Nz + j*Nz + k].ex=(  halfgrid_now[ (i+1)*Ny*Nz + j*Nz + k].ex + halfgrid_now[ (i+1)*Ny*Nz + j*Nz + k+1].ex +  halfgrid_now[ (i+1)*Ny*Nz + (j+2)*Nz + k].ex +  halfgrid_now[ (i+1)*Ny*Nz + (j+2)*Nz + k+2].ex)/4;
	//				grid_result[i*Ny*Nz + j*Nz + k].ey=(  halfgrid_now[ i*Ny*Nz + (j+1)*Nz + k].ey +  halfgrid_now[ (i+2)*Ny*Nz +( j+1) *Nz + k].ey + halfgrid_now[ i*Ny*Nz +( j+1) *Nz + k+2].ey+  halfgrid_now[ (i+2)*Ny*Nz +( j+1) *Nz + k+2].ey )/4;
	//				grid_result[i*Ny*Nz + j*Nz + k].ez= ( halfgrid_now[ i*Ny*Nz + j*Nz + k+1].ez + halfgrid_now[ i*Ny*Nz + (j+2)*Nz + k+1].ez + halfgrid_now[ (i+2)*Ny*Nz + j*Nz + k+1].ez + halfgrid_now[ (i+2)*Ny*Nz + (j+2)*Nz + k+1].ez)/4;
	//				grid_result[i*Ny*Nz + j*Nz + k].bx=(  halfgrid_now[ i*Ny*Nz + (j+1)*Nz + k+1].bx + halfgrid_now[ (i+2)*Ny*Nz + (j+1)*Nz + k+2].bx  )/2;
	//				grid_result[i*Ny*Nz + j*Nz + k].by=(  halfgrid_now[ (i+1)*Ny*Nz + j*Nz + k+1].by +  halfgrid_now[ (i+1)*Ny*Nz + (j+2)*Nz + k+1].by  )/2;
	//				grid_result[i*Ny*Nz + j*Nz + k].bz=(  halfgrid_now[ (i+1)*Ny*Nz + (j+1)*Nz + k ].bz +  halfgrid_now[ (i+1)*Ny*Nz + (j+1)*Nz + k+2].bz)/2;
	//			 }//else
	//		}	
	int result_x = 100;
	int result_y = 100;
	int result_z = 50;
	file << step << '\t' << halfgrid_now[result_x * Ny*Nz + result_y * Nz + result_z].ex << '\t' << halfgrid_now[result_x * Ny*Nz + result_y * Nz + result_z].ey << '\t' << halfgrid_now[result_x * Ny*Nz + result_y * Nz + result_z].ez << '\t';
	file << halfgrid_now[result_x * Ny*Nz + result_y * Nz + result_z].bx << '\t' << halfgrid_now[result_x * Ny*Nz + result_y * Nz + result_z].by << '\t' << halfgrid_now[result_x * Ny*Nz + result_y * Nz + result_z].bz << '\t';
	file << '\n';
	//save_result(halfgrid_now, step);
	step++;
	cout<<"Step--- "<<step<<" ---has finished."<<endl;
 }//while

}//函数结尾

void adi_fdtd_leapforg_mur(Grid* halfgrid_beforeX2,Grid* halfgrid_before, Grid* halfgrid_now)
{
	/*system("mkdir result");*/
	ofstream file("result\\leapforg_ADI_FDTD_steam0.1_mur02.txt");//用于保存结果
	//*******计算TE10模******//注意边界条件的问题,直接全部使用二阶mur吸收边界还是会出现问题，因为有棱边的存在，现在采用在棱边处使用一阶吸收边界条件，避免使用棱边

	int step = 0;//计算时间步长
	while (step< STEPS)
	{
		inject_field(halfgrid_before,halfgrid_now, step);
		//PART1---- 计算电场//
		for (int i = 0; i<Nx - 1; i++)
			for (int j = 1; j<Ny - 2; j++)
				for (int k = 1; k<Nz - 2; k++)
				{
			if (i == 0)//左边界，信号源处,一阶mur吸收边界处理棱边			
			{
				/*if (j == 0)
				{
					halfgrid_now[i*Ny*Nz + j*Nz + k].ex = 0;
					halfgrid_now[i*Ny*Nz + j*Nz + k].ey = 0;
					halfgrid_now[i*Ny*Nz + j*Nz + k].ez = halfgrid_before[(i+1) * Ny*Nz + (j+1) * Nz + k].ez + ((c*dt - sqrt(2)*dx) / (c*dt + sqrt(2)*dx))*(halfgrid_now[(i+1)* Ny*Nz + (j+1) * Nz + k].ez - halfgrid_before[i * Ny*Nz + j * Nz + k].ez);
				}
				else if (j == Ny - 1)
				{
					halfgrid_now[i*Ny*Nz + j*Nz + k].ex = 0;
					halfgrid_now[i*Ny*Nz + j*Nz + k].ey = 0;
					halfgrid_now[i*Ny*Nz + j*Nz + k].ez = halfgrid_before[(i+1) * Ny*Nz + (j-1) * Nz + k].ez + ((c*dt - sqrt(2)*dx) / (c*dt + sqrt(2)*dx))*(halfgrid_now[(i+1) * Ny*Nz + (j-1) * Nz + k].ez - halfgrid_before[i * Ny*Nz + j * Nz + k].ez);
				}
				else if (k == 0)
				{
					halfgrid_now[i*Ny*Nz + j*Nz + k].ex = 0;
					halfgrid_now[i*Ny*Nz + j*Nz + k].ez = 0;
					halfgrid_now[i*Ny*Nz + j*Nz + k].ey = halfgrid_before[(i + 1) * Ny*Nz + j * Nz + k+1].ey + ((c*dt - sqrt(2)*dx) / (c*dt + sqrt(2)*dx))*(halfgrid_now[(i + 1) * Ny*Nz + j * Nz + k+1].ey - halfgrid_before[i * Ny*Nz + j * Nz + k].ey);
				}
				else if (k == Nz - 1)
				{
					halfgrid_now[i*Ny*Nz + j*Nz + k].ex = 0;
					halfgrid_now[i*Ny*Nz + j*Nz + k].ez = 0;
					halfgrid_now[i*Ny*Nz + j*Nz + k].ey = halfgrid_before[(i + 1) * Ny*Nz + j * Nz + k - 1].ey + ((c*dt - sqrt(2)*dx) / (c*dt + sqrt(2)*dx))*(halfgrid_now[(i + 1) * Ny*Nz + j * Nz + k - 1].ey - halfgrid_before[i * Ny*Nz + j * Nz + k].ey);
				}*/
			  //else
				//{
					halfgrid_now[i*Ny*Nz + j*Nz + k].ex = 0;

					halfgrid_now[i*Ny*Nz + j*Nz + k].ey = -halfgrid_beforeX2[(i + 1)*Ny*Nz + j*Nz + k].ey + (c*dt - dx) / (c*dt + dx)*(halfgrid_now[(i + 1)*Ny*Nz + j*Nz + k].ey + halfgrid_beforeX2[i*Ny*Nz + j*Nz + k].ey)
					+ (2 * dx / (c*dt + dx))*(halfgrid_before[i*Ny*Nz + j*Nz + k].ey + halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k].ey)
					+ (dx*(c*dt)*(c*dt)) / (2 * dy*dy*(c*dt + dx))*(halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k].ey - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k].ey + halfgrid_before[i*Ny*Nz + (j - 1)*Nz + k].ey
					+ halfgrid_before[(i + 1)*Ny*Nz + (j + 1)*Nz + k].ey - 2 * halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k].ey + halfgrid_before[(i + 1)*Ny*Nz + (j - 1)*Nz + k].ey)
					+ (dx*(c*dt)*(c*dt)) / (2 * dz*dz*(c*dt + dx))*(halfgrid_before[i*Ny*Nz + j*Nz + k + 1].ey - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k].ey
					+ halfgrid_before[i*Ny*Nz + j*Nz + k - 1].ey - halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k + 1].ey - 2 * halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k].ey + halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k - 1].ey);
					
					halfgrid_now[i*Ny*Nz + j*Nz + k].ez = -halfgrid_beforeX2[(i + 1)*Ny*Nz + j*Nz + k].ez + (c*dt - dx) / (c*dt + dx)*(halfgrid_now[(i + 1)*Ny*Nz + j*Nz + k].ez + halfgrid_beforeX2[i*Ny*Nz + j*Nz + k].ez)
					+ (2 * dx / (c*dt + dx))*(halfgrid_before[i*Ny*Nz + j*Nz + k].ez + halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k].ez)
					+ (dx*(c*dt)*(c*dt)) / (2 * dy*dy*(c*dt + dx))*(halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k].ez - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k].ez + halfgrid_before[i*Ny*Nz + (j - 1)*Nz + k].ez
					+ halfgrid_before[(i + 1)*Ny*Nz + (j + 1)*Nz + k].ez - 2 * halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k].ez + halfgrid_before[(i + 1)*Ny*Nz + (j - 1)*Nz + k - 1].ez)
					+ (dx*(c*dt)*(c*dt)) / (2 * dz*dz*(c*dt + dx))*(halfgrid_before[i*Ny*Nz + j*Nz + k + 1].ez - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k].ez
					+ halfgrid_before[i*Ny*Nz + j*Nz + k - 1].ez - halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k + 1].ez - 2 * halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k].ez + halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k - 1].ez);
				//}

			}
			else if (i == Nx - 1)//右边界，信号源处
			{
					/*if (j == 0)
					{
						halfgrid_now[i*Ny*Nz + j*Nz + k].ex = 0;
						halfgrid_now[i*Ny*Nz + j*Nz + k].ey = 0;
						halfgrid_now[i*Ny*Nz + j*Nz + k].ez = halfgrid_before[(i - 1) * Ny*Nz + (j + 1) * Nz + k].ez + ((c*dt - sqrt(2)*dx) / (c*dt + sqrt(2)*dx))*(halfgrid_now[(i - 1)* Ny*Nz + (j + 1) * Nz + k].ez - halfgrid_before[i * Ny*Nz + j * Nz + k].ez);
					}
					else if (j == Ny - 1)
					{
						halfgrid_now[i*Ny*Nz + j*Nz + k].ex = 0;
						halfgrid_now[i*Ny*Nz + j*Nz + k].ey = 0;
						halfgrid_now[i*Ny*Nz + j*Nz + k].ez = halfgrid_before[(i - 1) * Ny*Nz + (j - 1) * Nz + k].ez + ((c*dt - sqrt(2)*dx) / (c*dt + sqrt(2)*dx))*(halfgrid_now[(i - 1) * Ny*Nz + (j - 1) * Nz + k].ez - halfgrid_before[i * Ny*Nz + j * Nz + k].ez);
					}
					else if (k == 0)
					{
						halfgrid_now[i*Ny*Nz + j*Nz + k].ex = 0;
						halfgrid_now[i*Ny*Nz + j*Nz + k].ez = 0;
						halfgrid_now[i*Ny*Nz + j*Nz + k].ey = halfgrid_before[(i - 1) * Ny*Nz + j * Nz + k + 1].ey + ((c*dt - sqrt(2)*dx) / (c*dt + sqrt(2)*dx))*(halfgrid_now[(i - 1) * Ny*Nz + j * Nz + k + 1].ey - halfgrid_before[i * Ny*Nz + j * Nz + k].ey);
					}
					else if (k == Nz - 1)
					{
						halfgrid_now[i*Ny*Nz + j*Nz + k].ex = 0;
						halfgrid_now[i*Ny*Nz + j*Nz + k].ez = 0;
						halfgrid_now[i*Ny*Nz + j*Nz + k].ey = halfgrid_before[(i - 1) * Ny*Nz + j * Nz + k - 1].ey + ((c*dt - sqrt(2)*dx) / (c*dt + sqrt(2)*dx))*(halfgrid_now[(i - 1) * Ny*Nz + j * Nz + k - 1].ey - halfgrid_before[i * Ny*Nz + j * Nz + k].ey);
					}
					else
					{*/
						halfgrid_now[i*Ny*Nz + j*Nz + k].ex = 0;
						halfgrid_now[i*Ny*Nz + j*Nz + k].ey = -halfgrid_beforeX2[(i - 1)*Ny*Nz + j*Nz + k].ey + (c*dt - dx) / (c*dt + dx)*(halfgrid_now[(i - 1)*Ny*Nz + j*Nz + k].ey + halfgrid_beforeX2[i*Ny*Nz + j*Nz + k].ey)
						+ (2 * dx / (c*dt + dx))*(halfgrid_before[i*Ny*Nz + j*Nz + k].ey + halfgrid_before[(i - 1)*Ny*Nz + j*Nz + k].ey)
						+ (dx*(c*dt)*(c*dt)) / (2 * dy*dy*(c*dt + dx))*(halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k].ey - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k].ey + halfgrid_before[i*Ny*Nz + (j - 1)*Nz + k].ey
					    + halfgrid_before[(i - 1)*Ny*Nz + j*Nz + k].ey - 2 * halfgrid_before[(i - 1)*Ny*Nz + j*Nz + k].ey + halfgrid_before[(i - 1)*Ny*Nz + (j - 1)*Nz + k].ey)
					    + (dx*(c*dt)*(c*dt)) / (2 * dz*dz*(c*dt + dx))*(halfgrid_before[i*Ny*Nz + j*Nz + k + 1].ey - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k].ey
					    + halfgrid_before[i*Ny*Nz + j*Nz + k - 1].ey - halfgrid_before[(i - 1)*Ny*Nz + j*Nz + k + 1].ey - 2 * halfgrid_before[(i - 1)*Ny*Nz + j*Nz + k].ey + halfgrid_before[(i - 1)*Ny*Nz + j*Nz + k - 1].ey);

				        halfgrid_now[i*Ny*Nz + j*Nz + k].ez = -halfgrid_beforeX2[(i - 1)*Ny*Nz + j*Nz + k].ez + (c*dt - dx) / (c*dt + dx)*(halfgrid_now[(i - 1)*Ny*Nz + j*Nz + k].ez + halfgrid_beforeX2[i*Ny*Nz + j*Nz + k].ez)
				    	+ (2 * dx / (c*dt + dx))*(halfgrid_before[i*Ny*Nz + j*Nz + k].ez + halfgrid_before[(i - 1)*Ny*Nz + j*Nz + k].ez)
					    + (dx*(c*dt)*(c*dt)) / (2 * dy*dy*(c*dt + dx))*(halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k].ez - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k].ez + halfgrid_before[i*Ny*Nz + (j - 1)*Nz + k].ez
					    + halfgrid_before[(i - 1)*Ny*Nz + j*Nz + k].ez - 2 * halfgrid_before[(i - 1)*Ny*Nz + j*Nz + k].ez + halfgrid_before[(i - 1)*Ny*Nz + (j - 1)*Nz + k].ez)
					    + (dx*(c*dt)*(c*dt)) / (2 * dz*dz*(c*dt + dx))*(halfgrid_before[i*Ny*Nz + j*Nz + k + 1].ez - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k].ez
					    + halfgrid_before[i*Ny*Nz + j*Nz + k - 1].ez - halfgrid_before[(i - 1)*Ny*Nz + j*Nz + k + 1].ez - 2 * halfgrid_before[(i - 1)*Ny*Nz + j*Nz + k].ez + halfgrid_before[(i - 1)*Ny*Nz + j*Nz + k - 1].ez);

					//}
															
			}
			/*else if (j == 0)
			{
				if (i == 0)
				{
					halfgrid_now[i*Ny*Nz + j*Nz + k].ex = 0;
					halfgrid_now[i*Ny*Nz + j*Nz + k].ey = 0;
					halfgrid_now[i*Ny*Nz + j*Nz + k].ez = halfgrid_before[(i + 1) * Ny*Nz + (j + 1) * Nz + k].ez + ((c*dt - sqrt(2)*dy) / (c*dt + sqrt(2)*dy))*(halfgrid_now[(i + 1)* Ny*Nz + (j + 1) * Nz + k].ez - halfgrid_before[i * Ny*Nz + j * Nz + k].ez);
				}
				else if (i == Nx - 1)
				{
					halfgrid_now[i*Ny*Nz + j*Nz + k].ex = 0;
					halfgrid_now[i*Ny*Nz + j*Nz + k].ey = 0;
					halfgrid_now[i*Ny*Nz + j*Nz + k].ez = halfgrid_before[(i - 1) * Ny*Nz + (j + 1) * Nz + k].ez + ((c*dt - sqrt(2)*dy) / (c*dt + sqrt(2)*dy))*(halfgrid_now[(i - 1) * Ny*Nz + (j + 1) * Nz + k].ez - halfgrid_before[i * Ny*Nz + j * Nz + k].ez);
				}
				else if (k == 0)
				{
					halfgrid_now[i*Ny*Nz + j*Nz + k].ey = 0;
					halfgrid_now[i*Ny*Nz + j*Nz + k].ez = 0;
					halfgrid_now[i*Ny*Nz + j*Nz + k].ex = halfgrid_before[i * Ny*Nz + (j + 1) * Nz + k + 1].ex + ((c*dt - sqrt(2)*dy) / (c*dt + sqrt(2)*dy))*(halfgrid_now[(i - 1) * Ny*Nz + (j + 1) * Nz + k + 1].ex - halfgrid_before[i * Ny*Nz + j * Nz + k].ex);
				}
				else if (k == Nz - 1)
				{
					halfgrid_now[i*Ny*Nz + j*Nz + k].ey = 0;
					halfgrid_now[i*Ny*Nz + j*Nz + k].ez = 0;
					halfgrid_now[i*Ny*Nz + j*Nz + k].ex = halfgrid_before[i * Ny*Nz + (j+1) * Nz + k - 1].ex + ((c*dt - sqrt(2)*dy) / (c*dt + sqrt(2)*dy))*(halfgrid_now[i * Ny*Nz + (j+1) * Nz + k - 1].ex - halfgrid_before[i * Ny*Nz + j * Nz + k].ex);
				}
				else
				{				
			    halfgrid_now[i*Ny*Nz + j*Nz + k].ey = 0;

				halfgrid_now[i*Ny*Nz + j*Nz + k].ex = -halfgrid_beforeX2[i*Ny*Nz + (j + 1)*Nz + k].ex + (c*dt - dy) / (c*dt + dy)*(halfgrid_now[i*Ny*Nz + (j + 1)*Nz + k].ex + halfgrid_beforeX2[i*Ny*Nz + j*Nz + k].ex)
					+ (2 * dy / (c*dt + dy))*(halfgrid_before[i*Ny*Nz + j*Nz + k].ex + halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k].ex)
					+ (dy*(c*dt)*(c*dt)) / (2 * dz*dz*(c*dt + dy))*(halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k].ex - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k].ex + halfgrid_before[(i - 1)*Ny*Nz + j*Nz + k].ex
					+ halfgrid_before[(i + 1)*Ny*Nz + (j + 1)*Nz + k].ex - 2 * halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k].ex + halfgrid_before[(i - 1)*Ny*Nz + (j + 1)*Nz + k].ex)
					+ (dx*(c*dt)*(c*dt)) / (2 * dz*dz*(c*dt + dx))*(halfgrid_before[i*Ny*Nz + j*Nz + k + 1].ex - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k].ex
					+ halfgrid_before[i*Ny*Nz + j*Nz + k - 1].ex - halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k + 1].ex - 2 * halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k].ex + halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k - 1].ex);

				halfgrid_now[i*Ny*Nz + j*Nz + k].ez = -halfgrid_beforeX2[i*Ny*Nz + (j + 1)*Nz + k].ez + (c*dt - dy) / (c*dt + dy)*(halfgrid_now[i*Ny*Nz + (j + 1)*Nz + k].ez + halfgrid_beforeX2[i*Ny*Nz + j*Nz + k].ez)
					+ (2 * dy / (c*dt + dy))*(halfgrid_before[i*Ny*Nz + j*Nz + k].ez + halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k].ez)
					+ (dy*(c*dt)*(c*dt)) / (2 * dz*dz*(c*dt + dy))*(halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k].ez - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k].ez + halfgrid_before[(i - 1)*Ny*Nz + j*Nz + k].ez
					+ halfgrid_before[(i + 1)*Ny*Nz + (j + 1)*Nz + k].ez - 2 * halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k].ez + halfgrid_before[(i - 1)*Ny*Nz + (j + 1)*Nz + k].ez)
					+ (dx*(c*dt)*(c*dt)) / (2 * dz*dz*(c*dt + dx))*(halfgrid_before[i*Ny*Nz + j*Nz + k + 1].ez - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k].ez
					+ halfgrid_before[i*Ny*Nz + j*Nz + k - 1].ez - halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k + 1].ez - 2 * halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k].ez + halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k - 1].ez);

				}
				
			}
			else if (j == Ny - 1)
			{
				if (i == 0)
				{
					halfgrid_now[i*Ny*Nz + j*Nz + k].ex = 0;
					halfgrid_now[i*Ny*Nz + j*Nz + k].ey = 0;
					halfgrid_now[i*Ny*Nz + j*Nz + k].ez = halfgrid_before[(i + 1) * Ny*Nz + (j - 1) * Nz + k].ez + ((c*dt - sqrt(2)*dy) / (c*dt + sqrt(2)*dy))*(halfgrid_now[(i + 1)* Ny*Nz + (j - 1) * Nz + k].ez - halfgrid_before[i * Ny*Nz + j * Nz + k].ez);
				}
				else if (i == Nx - 1)
				{
					halfgrid_now[i*Ny*Nz + j*Nz + k].ex = 0;
					halfgrid_now[i*Ny*Nz + j*Nz + k].ey = 0;
					halfgrid_now[i*Ny*Nz + j*Nz + k].ez = halfgrid_before[(i - 1) * Ny*Nz + (j - 1) * Nz + k].ez + ((c*dt - sqrt(2)*dy) / (c*dt + sqrt(2)*dy))*(halfgrid_now[(i - 1) * Ny*Nz + (j - 1) * Nz + k].ez - halfgrid_before[i * Ny*Nz + j * Nz + k].ez);
				}
				else if (k == 0)
				{
					halfgrid_now[i*Ny*Nz + j*Nz + k].ey = 0;
					halfgrid_now[i*Ny*Nz + j*Nz + k].ez = 0;
					halfgrid_now[i*Ny*Nz + j*Nz + k].ex = halfgrid_before[i * Ny*Nz + (j - 1) * Nz + k + 1].ex + ((c*dt - sqrt(2)*dy) / (c*dt + sqrt(2)*dy))*(halfgrid_now[(i - 1) * Ny*Nz + (j - 1) * Nz + k + 1].ex - halfgrid_before[i * Ny*Nz + j * Nz + k].ex);
				}
				else if (k == Nz - 1)
				{
					halfgrid_now[i*Ny*Nz + j*Nz + k].ey = 0;
					halfgrid_now[i*Ny*Nz + j*Nz + k].ez = 0;
					halfgrid_now[i*Ny*Nz + j*Nz + k].ex = halfgrid_before[i * Ny*Nz + (j - 1) * Nz + k - 1].ex + ((c*dt - sqrt(2)*dy) / (c*dt + sqrt(2)*dy))*(halfgrid_now[i * Ny*Nz + (j - 1) * Nz + k - 1].ex - halfgrid_before[i * Ny*Nz + j * Nz + k].ex);
				}
				else
				{
				    halfgrid_now[i*Ny*Nz + j*Nz + k].ey = 0;

				    halfgrid_now[i*Ny*Nz + j*Nz + k].ex = -halfgrid_beforeX2[i*Ny*Nz + (j - 1)*Nz + k].ex + (c*dt - dy) / (c*dt + dy)*(halfgrid_now[i*Ny*Nz + (j - 1)*Nz + k].ex + halfgrid_beforeX2[i*Ny*Nz + j*Nz + k].ex)
					+ (2 * dy / (c*dt + dy))*(halfgrid_before[i*Ny*Nz + j*Nz + k].ex + halfgrid_before[i*Ny*Nz + (j - 1)*Nz + k].ex)
					+ (dy*(c*dt)*(c*dt)) / (2 * dz*dz*(c*dt + dy))*(halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k].ex - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k].ex + halfgrid_before[(i - 1)*Ny*Nz + j*Nz + k].ex
					+ halfgrid_before[(i + 1)*Ny*Nz + (j - 1)*Nz + k].ex - 2 * halfgrid_before[i*Ny*Nz + (j - 1)*Nz + k].ex + halfgrid_before[(i - 1)*Ny*Nz + (j - 1)*Nz + k].ex)
					+ (dx*(c*dt)*(c*dt)) / (2 * dz*dz*(c*dt + dx))*(halfgrid_before[i*Ny*Nz + j*Nz + k + 1].ex - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k].ex
					+ halfgrid_before[i*Ny*Nz + j*Nz + k - 1].ex - halfgrid_before[i*Ny*Nz + (j - 1)*Nz + k + 1].ex - 2 * halfgrid_before[i*Ny*Nz + (j - 1)*Nz + k].ex + halfgrid_before[i*Ny*Nz + (j - 1)*Nz + k - 1].ex);

				    halfgrid_now[i*Ny*Nz + j*Nz + k].ez = -halfgrid_beforeX2[i*Ny*Nz + (j - 1)*Nz + k].ez + (c*dt - dy) / (c*dt + dy)*(halfgrid_now[i*Ny*Nz + (j - 1)*Nz + k].ez + halfgrid_beforeX2[i*Ny*Nz + j*Nz + k].ez)
					+ (2 * dy / (c*dt + dy))*(halfgrid_before[i*Ny*Nz + j*Nz + k].ez + halfgrid_before[i*Ny*Nz + (j - 1)*Nz + k].ez)
					+ (dy*(c*dt)*(c*dt)) / (2 * dz*dz*(c*dt + dy))*(halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k].ez - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k].ez + halfgrid_before[(i - 1)*Ny*Nz + j*Nz + k].ez
					+ halfgrid_before[(i + 1)*Ny*Nz + (j - 1)*Nz + k].ez - 2 * halfgrid_before[i*Ny*Nz + (j - 1)*Nz + k].ez + halfgrid_before[(i - 1)*Ny*Nz + (j - 1)*Nz + k].ez)
					+ (dx*(c*dt)*(c*dt)) / (2 * dz*dz*(c*dt + dx))*(halfgrid_before[i*Ny*Nz + j*Nz + k + 1].ez - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k].ez
					+ halfgrid_before[i*Ny*Nz + j*Nz + k - 1].ez - halfgrid_before[i*Ny*Nz + (j - 1)*Nz + k + 1].ez - 2 * halfgrid_before[i*Ny*Nz + (j - 1)*Nz + k].ez + halfgrid_before[i*Ny*Nz + (j - 1)*Nz + k - 1].ez);

				}
			}

			else if (k == 0)
			{
				if (i == 0)
				{
					halfgrid_now[i*Ny*Nz + j*Nz + k].ex = 0;
					halfgrid_now[i*Ny*Nz + j*Nz + k].ez = 0;
					halfgrid_now[i*Ny*Nz + j*Nz + k].ey = halfgrid_before[(i + 1) * Ny*Nz + j * Nz + k + 1].ey + ((c*dt - sqrt(2)*dz) / (c*dt + sqrt(2)*dz))*(halfgrid_now[(i + 1)* Ny*Nz + j * Nz + k + 1].ey - halfgrid_before[i * Ny*Nz + j * Nz + k].ey);
				}
				else if (i == Nx - 1)
				{
					halfgrid_now[i*Ny*Nz + j*Nz + k].ex = 0;
					halfgrid_now[i*Ny*Nz + j*Nz + k].ez = 0;
					halfgrid_now[i*Ny*Nz + j*Nz + k].ey = halfgrid_before[(i - 1) * Ny*Nz + j * Nz + k + 1].ey + ((c*dt - sqrt(2)*dz) / (c*dt + sqrt(2)*dz))*(halfgrid_now[(i - 1) * Ny*Nz + j * Nz + k+1].ey - halfgrid_before[i * Ny*Nz + j * Nz + k].ey);
				}
				else if (j == 0)
				{
					halfgrid_now[i*Ny*Nz + j*Nz + k].ey = 0;
					halfgrid_now[i*Ny*Nz + j*Nz + k].ez = 0;
					halfgrid_now[i*Ny*Nz + j*Nz + k].ex = halfgrid_before[i * Ny*Nz + (j + 1) * Nz + k + 1].ex + ((c*dt - sqrt(2)*dz) / (c*dt + sqrt(2)*dz))*(halfgrid_now[i * Ny*Nz + (j + 1) * Nz + k + 1].ex - halfgrid_before[i * Ny*Nz + j * Nz + k].ex);
				}
				else if (j == Ny - 1)
				{
					halfgrid_now[i*Ny*Nz + j*Nz + k].ey = 0;
					halfgrid_now[i*Ny*Nz + j*Nz + k].ez = 0;
					halfgrid_now[i*Ny*Nz + j*Nz + k].ex = halfgrid_before[i * Ny*Nz + (j - 1) * Nz + k + 1].ex + ((c*dt - sqrt(2)*dz) / (c*dt + sqrt(2)*dz))*(halfgrid_now[i * Ny*Nz + (j - 1) * Nz + k + 1].ex - halfgrid_before[i * Ny*Nz + j * Nz + k].ex);
				}
				else
				{
					halfgrid_now[i*Ny*Nz + j*Nz + k].ez = 0;

				    halfgrid_now[i*Ny*Nz + j*Nz + k].ex = -halfgrid_beforeX2[i*Ny*Nz + j*Nz + k + 1].ex + (c*dt - dy) / (c*dt + dy)*(halfgrid_now[i*Ny*Nz + j*Nz + k + 1].ex + halfgrid_beforeX2[i*Ny*Nz + j*Nz + k].ex)
					+ (2 * dz / (c*dt + dz))*(halfgrid_before[i*Ny*Nz + j*Nz + k].ex + halfgrid_before[i*Ny*Nz + j*Nz + k + 1].ex)
					+ (dz*(c*dt)*(c*dt)) / (2 * dx*dx*(c*dt + dz))*(halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k].ex - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k].ex + halfgrid_before[(i - 1)*Ny*Nz + j*Nz + k].ex
					+ halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k + 1].ex - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k + 1].ex + halfgrid_before[i*Ny*Nz + j*Nz + k + 1].ex)
					+ (dz*(c*dt)*(c*dt)) / (2 * dy*dy*(c*dt + dz))*(halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k].ex - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k].ex
					+ halfgrid_before[i*Ny*Nz + (j - 1)*Nz + k].ex - halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k + 1].ex - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k + 1].ex + halfgrid_before[i*Ny*Nz + (j - 1)*Nz + k + 1].ex);

				    halfgrid_now[i*Ny*Nz + j*Nz + k].ey = -halfgrid_beforeX2[i*Ny*Nz + j*Nz + k + 1].ey + (c*dt - dy) / (c*dt + dy)*(halfgrid_now[i*Ny*Nz + j*Nz + k + 1].ey + halfgrid_beforeX2[i*Ny*Nz + j*Nz + k].ey)
					+ (2 * dz / (c*dt + dz))*(halfgrid_before[i*Ny*Nz + j*Nz + k].ey + halfgrid_before[i*Ny*Nz + j*Nz + k + 1].ey)
					+ (dz*(c*dt)*(c*dt)) / (2 * dx*dx*(c*dt + dz))*(halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k].ey - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k].ey + halfgrid_before[(i - 1)*Ny*Nz + j*Nz + k].ey
					+ halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k + 1].ey - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k + 1].ey + halfgrid_before[i*Ny*Nz + j*Nz + k + 1].ey)
					+ (dz*(c*dt)*(c*dt)) / (2 * dy*dy*(c*dt + dz))*(halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k].ey - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k].ey
					+ halfgrid_before[i*Ny*Nz + (j - 1)*Nz + k].ey - halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k + 1].ey - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k + 1].ey + halfgrid_before[i*Ny*Nz + (j - 1)*Nz + k + 1].ey);
				}

				

			}
			else if (k == Nz - 1)
			{
				if (i == 0)
				{
					halfgrid_now[i*Ny*Nz + j*Nz + k].ex = 0;
					halfgrid_now[i*Ny*Nz + j*Nz + k].ez = 0;
					halfgrid_now[i*Ny*Nz + j*Nz + k].ey = halfgrid_before[(i + 1) * Ny*Nz + j * Nz + k - 1].ey + ((c*dt - sqrt(2)*dz) / (c*dt + sqrt(2)*dz))*(halfgrid_now[(i + 1)* Ny*Nz + j * Nz + k - 1].ey - halfgrid_before[i * Ny*Nz + j * Nz + k].ey);
				}
				else if (i == Nx - 1)
				{
					halfgrid_now[i*Ny*Nz + j*Nz + k].ex = 0;
					halfgrid_now[i*Ny*Nz + j*Nz + k].ez = 0;
					halfgrid_now[i*Ny*Nz + j*Nz + k].ey = halfgrid_before[(i - 1) * Ny*Nz + j * Nz + k - 1].ey + ((c*dt - sqrt(2)*dz) / (c*dt + sqrt(2)*dz))*(halfgrid_now[(i - 1) * Ny*Nz + j * Nz + k - 1].ey - halfgrid_before[i * Ny*Nz + j * Nz + k].ey);
				}
				else if (j == 0)
				{
					halfgrid_now[i*Ny*Nz + j*Nz + k].ey = 0;
					halfgrid_now[i*Ny*Nz + j*Nz + k].ez = 0;
					halfgrid_now[i*Ny*Nz + j*Nz + k].ex = halfgrid_before[i * Ny*Nz + (j + 1) * Nz + k - 1].ex + ((c*dt - sqrt(2)*dz) / (c*dt + sqrt(2)*dz))*(halfgrid_now[i * Ny*Nz + (j + 1) * Nz + k - 1].ex - halfgrid_before[i * Ny*Nz + j * Nz + k].ex);
				}
				else if (j == Ny - 1)
				{
					halfgrid_now[i*Ny*Nz + j*Nz + k].ey = 0;
					halfgrid_now[i*Ny*Nz + j*Nz + k].ez = 0;
					halfgrid_now[i*Ny*Nz + j*Nz + k].ex = halfgrid_before[i * Ny*Nz + (j - 1) * Nz + k - 1].ex + ((c*dt - sqrt(2)*dz) / (c*dt + sqrt(2)*dz))*(halfgrid_now[i * Ny*Nz + (j - 1) * Nz + k - 1].ex - halfgrid_before[i * Ny*Nz + j * Nz + k].ex);
				}
				else
				{
					halfgrid_now[i*Ny*Nz + j*Nz + k].ez = 0;

				    halfgrid_now[i*Ny*Nz + j*Nz + k].ex = -halfgrid_beforeX2[i*Ny*Nz + j*Nz + k - 1].ex + (c*dt - dy) / (c*dt + dy)*(halfgrid_now[i*Ny*Nz + j*Nz + k - 1].ex + halfgrid_beforeX2[i*Ny*Nz + j*Nz + k].ex)
					+ (2 * dz / (c*dt + dz))*(halfgrid_before[i*Ny*Nz + j*Nz + k].ex + halfgrid_before[i*Ny*Nz + j*Nz + k - 1].ex)
					+ (dz*(c*dt)*(c*dt)) / (2 * dx*dx*(c*dt + dz))*(halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k].ex - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k].ex + halfgrid_before[(i - 1)*Ny*Nz + j*Nz + k].ex
					+ halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k - 1].ex - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k - 1].ex + halfgrid_before[i*Ny*Nz + j*Nz + k - 1].ex)
					+ (dz*(c*dt)*(c*dt)) / (2 * dy*dy*(c*dt + dz))*(halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k].ex - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k].ex
					+ halfgrid_before[i*Ny*Nz + (j - 1)*Nz + k].ex - halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k - 1].ex - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k - 1].ex + halfgrid_before[i*Ny*Nz + (j - 1)*Nz + k - 1].ex);

				    halfgrid_now[i*Ny*Nz + j*Nz + k].ey = -halfgrid_beforeX2[i*Ny*Nz + j*Nz + k - 1].ey + (c*dt - dy) / (c*dt + dy)*(halfgrid_now[i*Ny*Nz + j*Nz + k - 1].ey + halfgrid_beforeX2[i*Ny*Nz + j*Nz + k].ey)
					+ (2 * dz / (c*dt + dz))*(halfgrid_before[i*Ny*Nz + j*Nz + k].ey + halfgrid_before[i*Ny*Nz + j*Nz + k - 1].ey)
					+ (dz*(c*dt)*(c*dt)) / (2 * dx*dx*(c*dt + dz))*(halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k].ey - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k].ey + halfgrid_before[(i - 1)*Ny*Nz + j*Nz + k].ey
					+ halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k - 1].ey - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k - 1].ey + halfgrid_before[i*Ny*Nz + j*Nz + k - 1].ey)
					+ (dz*(c*dt)*(c*dt)) / (2 * dy*dy*(c*dt + dz))*(halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k].ey - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k].ey
					+ halfgrid_before[i*Ny*Nz + (j - 1)*Nz + k].ey - halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k - 1].ey - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k - 1].ey + halfgrid_before[i*Ny*Nz + (j - 1)*Nz + k - 1].ey);

				}
				
			}
			*/

			else
			{

				//ex
				double ty0 = (halfgrid_before[i*Ny*Nz + j*Nz + k].bz - halfgrid_before[i*Ny*Nz + (j - 1)*Nz + k].bz) / dy;
				double tz0 = (halfgrid_before[i*Ny*Nz + j*Nz + k].by - halfgrid_before[i*Ny*Nz + j*Nz + k - 1].by) / dz;
				halfgrid_now[i*Ny*Nz + j*Nz + k].ex = halfgrid_before[i*Ny*Nz + j*Nz + k].ex + dt*(1 / epsl_x)*(ty0 - tz0);
				//ey
				double tx1 = (halfgrid_before[i*Ny*Nz + j*Nz + k].bx - halfgrid_before[i*Ny*Nz + j*Nz + k - 1].bx) / dz;
				double tz1 = (halfgrid_before[i*Ny*Nz + j*Nz + k].bz - halfgrid_before[(i - 1)*Ny*Nz + j*Nz + k].bz) / dx;
				halfgrid_now[i*Ny*Nz + j*Nz + k].ey = halfgrid_before[i*Ny*Nz + j*Nz + k].ey + dt*(1 / epsl_y)*(tx1 - tz1);
				//ez
				double ty2 = (halfgrid_before[i*Ny*Nz + j*Nz + k].by - halfgrid_before[(i - 1)*Ny*Nz + j*Nz + k].by) / dx;
				double tx2 = (halfgrid_before[i*Ny*Nz + j*Nz + k].bx - halfgrid_before[i*Ny*Nz + (j - 1)*Nz + k].bx) / dy;
				halfgrid_now[i*Ny*Nz + j*Nz + k].ez = halfgrid_before[i*Ny*Nz + j*Nz + k].ez + dt*(1 / epsl_z)*(ty2 - tx2);

			}//else
			halfgrid_beforeX2[i*Ny*Nz + j*Nz + k].ex = halfgrid_before[i*Ny*Nz + j*Nz + k].ex;//保存上上步的结果
			halfgrid_beforeX2[i*Ny*Nz + j*Nz + k].ey = halfgrid_before[i*Ny*Nz + j*Nz + k].ey;
			halfgrid_beforeX2[i*Ny*Nz + j*Nz + k].ez = halfgrid_before[i*Ny*Nz + j*Nz + k].ez;


			halfgrid_before[i*Ny*Nz + j*Nz + k].ex = halfgrid_now[i*Ny*Nz + j*Nz + k].ex;//保存前一步结果
			halfgrid_before[i*Ny*Nz + j*Nz + k].ey = halfgrid_now[i*Ny*Nz + j*Nz + k].ey;
			halfgrid_before[i*Ny*Nz + j*Nz + k].ez = halfgrid_now[i*Ny*Nz + j*Nz + k].ez;
				}

		//PART2----计算磁场//
		for (int i = 0; i<Nx - 1; i++)
			for (int j = 0; j<Ny - 1; j++)
				for (int k = 0; k < Nz - 1; k++)
				{

			//bx
			double ty3 = (halfgrid_before[i*Ny*Nz + j*Nz + k + 1].ey - halfgrid_before[i*Ny*Nz + j*Nz + k].ey) / dz;
			double tz3 = (halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k].ez - halfgrid_before[i*Ny*Nz + j*Nz + k].ez) / dy;
			halfgrid_now[i*Ny*Nz + j*Nz + k].bx = halfgrid_before[i*Ny*Nz + j*Nz + k].bx + dt*(1 / mur_x)*(ty3 - tz3);
			//by
			double tz4 = (halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k].ez - halfgrid_before[i*Ny*Nz + j*Nz + k].ez) / dx;
			double tx4 = (halfgrid_before[i*Ny*Nz + j*Nz + k + 1].ex - halfgrid_before[i*Ny*Nz + j*Nz + k].ex) / dz;
			halfgrid_now[i*Ny*Nz + j*Nz + k].by = halfgrid_before[i*Ny*Nz + j*Nz + k].by + dt*(1 / mur_y)*(tz4 - tx4);
			//bz
			double tx5 = (halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k].ex - halfgrid_before[i*Ny*Nz + j*Nz + k].ex) / dy;
			double ty5 = (halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k].ey - halfgrid_before[i*Ny*Nz + j*Nz + k].ey) / dx;
			halfgrid_now[i*Ny*Nz + j*Nz + k].bz = halfgrid_before[i*Ny*Nz + j*Nz + k].bz + dt*(1 / mur_z)*(tx5 - ty5);

			halfgrid_beforeX2[i*Ny*Nz + j*Nz + k].bx = halfgrid_before[i*Ny*Nz + j*Nz + k].bx;//保存上上步的结果
			halfgrid_beforeX2[i*Ny*Nz + j*Nz + k].by = halfgrid_before[i*Ny*Nz + j*Nz + k].by;
			halfgrid_beforeX2[i*Ny*Nz + j*Nz + k].bz = halfgrid_before[i*Ny*Nz + j*Nz + k].bz;

			halfgrid_before[i*Ny*Nz + j*Nz + k].bx = halfgrid_now[i*Ny*Nz + j*Nz + k].bx;//保存为前一步结果
			halfgrid_before[i*Ny*Nz + j*Nz + k].by = halfgrid_now[i*Ny*Nz + j*Nz + k].by;
			halfgrid_before[i*Ny*Nz + j*Nz + k].bz = halfgrid_now[i*Ny*Nz + j*Nz + k].bz;
				}

		////若取平均，则网格点处的电磁场计算如下
		//
		//for(int i = 0; i<Nx;i++)
		//	for(int j =0; j<Ny; j++)
		//		for (int k = 0; k < Nz; k++)
		//		{
		//			  if(i == 0 )//左边界，信号源处
		//			 {
		//				
		//			 }
		//		           
		//		     if(i== Nx - 1 )//右边界，吸收边界
		//			 {
		//				
		//			 }
		//			 if(k == Nz - 1|| k==0 || j==Ny - 1 || j==0 )//四个平面
		//			 {
		//				 
		//			 }
		//			 else
		//			 {
		//			 
		//			    grid_result[i*Ny*Nz + j*Nz + k].ex=(  halfgrid_now[ (i+1)*Ny*Nz + j*Nz + k].ex + halfgrid_now[ (i+1)*Ny*Nz + j*Nz + k+1].ex +  halfgrid_now[ (i+1)*Ny*Nz + (j+2)*Nz + k].ex +  halfgrid_now[ (i+1)*Ny*Nz + (j+2)*Nz + k+2].ex)/4;
		//				grid_result[i*Ny*Nz + j*Nz + k].ey=(  halfgrid_now[ i*Ny*Nz + (j+1)*Nz + k].ey +  halfgrid_now[ (i+2)*Ny*Nz +( j+1) *Nz + k].ey + halfgrid_now[ i*Ny*Nz +( j+1) *Nz + k+2].ey+  halfgrid_now[ (i+2)*Ny*Nz +( j+1) *Nz + k+2].ey )/4;
		//				grid_result[i*Ny*Nz + j*Nz + k].ez= ( halfgrid_now[ i*Ny*Nz + j*Nz + k+1].ez + halfgrid_now[ i*Ny*Nz + (j+2)*Nz + k+1].ez + halfgrid_now[ (i+2)*Ny*Nz + j*Nz + k+1].ez + halfgrid_now[ (i+2)*Ny*Nz + (j+2)*Nz + k+1].ez)/4;
		//				grid_result[i*Ny*Nz + j*Nz + k].bx=(  halfgrid_now[ i*Ny*Nz + (j+1)*Nz + k+1].bx + halfgrid_now[ (i+2)*Ny*Nz + (j+1)*Nz + k+2].bx  )/2;
		//				grid_result[i*Ny*Nz + j*Nz + k].by=(  halfgrid_now[ (i+1)*Ny*Nz + j*Nz + k+1].by +  halfgrid_now[ (i+1)*Ny*Nz + (j+2)*Nz + k+1].by  )/2;
		//				grid_result[i*Ny*Nz + j*Nz + k].bz=(  halfgrid_now[ (i+1)*Ny*Nz + (j+1)*Nz + k ].bz +  halfgrid_now[ (i+1)*Ny*Nz + (j+1)*Nz + k+2].bz)/2;
		//			 }//else
		//		}	
		int result_x = 100; 
		int result_y = 100;
		int result_z = 50;
		file << step << '\t' << halfgrid_now[result_x * Ny*Nz + result_y * Nz + result_z].ex << '\t' << halfgrid_now[result_x * Ny*Nz + result_y * Nz + result_z].ey << '\t' << halfgrid_now[result_x * Ny*Nz + result_y * Nz + result_z].ez << '\t';
		file << halfgrid_now[result_x * Ny*Nz + result_y * Nz + result_z].bx << '\t' << halfgrid_now[result_x * Ny*Nz + result_y * Nz + result_z].by << '\t' << halfgrid_now[result_x * Ny*Nz + result_y * Nz + result_z].bz << '\t';
		file << '\n';
		//save_result(halfgrid_now, step);
		step++;
		cout << "Step--- " << step << " ---has finished." << endl;
	}//while

}//函数结尾


#endif