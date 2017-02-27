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
		cout << "请输入吸收边界类型:\n" << endl;
		
		cout << "1. 金属吸收边界\n" << endl;
		cout << "2. Mur吸收边界\n" << endl;
		//cout << "3. 二阶Mur吸收边界\n" << endl;
		cin >> choice;
		switch (choice)
		{

		case 0:
		{
			fdtd_matel(halfgrid_now);//进行计算
			break;
		}
		case 1:
		{
			adi_fdtd_leapforg_matel_GSS(halfgrid_now);//进行计算
			break;
		}
		case 2:
		{
			adi_fdtd_leapforg_mur1_GSS(halfgrid_now);//进行计算
			break;
		}
		//case 3:
		//{
		//	adi_fdtd_leapforg_mur2_GSS(halfgrid_beforeX2, halfgrid_before, halfgrid_now);//进行计算
		//	break;
		//}
		default:
		{
			cout << "抱歉，无此选项，请重新输入。\n" << endl;
			continue;
		}

		}

		free(halfgrid_beforeX2, halfgrid_before, halfgrid_now);//保存结果,并释放空间

		clock_t end_time = clock();
		cout << "Running time: " << static_cast<float>(end_time - start_time) / CLOCKS_PER_SEC << " s" << endl;
		system("pause");
		return 0;
	}
}