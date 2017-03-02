//writen by liuzhichao 
//2016/10/8
#include "definer.h"
#include "grid.h"
#include "functions.h"
#include "GSS_ADI_matel.h"
#include "GSS_ADI_mur_1st.h"
#include "GSS_ADI_mur_2nd.h"
#include "3D_FDTD.h"


int main()
{
	clock_t start_time = clock();
	initGrid(halfgrid_beforeX2, halfgrid_before, halfgrid_now);//对网格数据进行初始化
	int choice;

	while (TRUE)
	{
		cout << ">>>>>>>>>>>>>请 选 择:<<<<<<<<<<<<\n" << endl;
		
		cout << ">1. 金 属 吸 收 边 界\n" << endl;
		cout << ">2. Mur 吸 收 边 界\n" << endl;
		cout << ">3. 计 算 CFL 稳 定 性 条 件 dt\n" << endl;
		cout << ">4. 计 算 理 论 值\n" << endl;
		cout << ">5. 结 束 程 序\n" << endl;
		cin >> choice;
		cout << '\n' << endl;
		switch (choice)
		{

		case 0:
		{
			fdtd_matel(halfgrid_now);//进行计算
			continue;
		}
		case 1:
		{
			adi_fdtd_leapforg_matel_GSS(halfgrid_now);//进行计算
			continue;
		}
		case 2:
		{
			adi_fdtd_leapforg_mur1_GSS(halfgrid_now);//进行计算
			continue;
		}
		case 3:
		{
			double result = CFL_calc();
			cout <<"在满足CFL稳定性条件下: dt <= "<< result << endl;
			cout << '\n' << endl;//进行计算
			continue;
		}
		case 4:
		{
			theor_val_gen();
			cout << "理论值已产生并存储至---> "<< theor_val_filepath << "\n" << endl;
			continue;

		}
		case 5:
		{
			cout << ">>>>>>>>程序结束<<<<<<<<" << endl;
			break;

		}
		default:
		{
			cout << "抱歉，无此选项，请重新输入。\n" << endl;
			continue;
		}

		}

		free(halfgrid_beforeX2, halfgrid_before, halfgrid_now);//释放空间

		clock_t end_time = clock();
		cout << "Running time: " << static_cast<float>(end_time - start_time) / CLOCKS_PER_SEC << " s" << endl;
		system("pause");
		return 0;
	}
}