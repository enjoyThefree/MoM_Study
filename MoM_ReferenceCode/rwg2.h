#pragma once
/***********************************************************/
/*            该函数是对九点积分法的数据预处理             */
/*                 得到计算中需要一些矢量                  */
/***********************************************************/
/*                                                         */
/*                   Author：Yanlin Xu                     */
/*                                                         */
/***********************************************************/
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "rwg2.h"
#include "armadillo"

using namespace arma;
using namespace std;

extern int TrianglesTotal, EdgesTotal;
extern mat t_data;
extern mat p_data;
extern mat EdgeElement;
extern vec TrianglePlus;
extern vec TriangleMinus;
extern vec EdgeLength;
extern mat Center;
extern vec Area;

mat  Pho_Plus;
mat  Pho_Minus;
cube Center_;     
cube Pho_Plus_;     //对于正element三角形 九点划分后 非公共边点到 小三角形重心的向量坐标
cube Pho_Minus_;

void rwg2_main()
{
	rwg1_main();

	int n1 = 0, n2 = 0, n3 = 0;
	int NoPlus = 0, NoMinus = 0, NODE = 0;

	rowvec r1(1, 3), r2(1, 3), r3(1, 3);
	rowvec r12(1, 3), r23(1, 3), r13(1, 3);
	rowvec M(1, 3), C1(1, 3), C2(1, 3), C3(1, 3), C4(1, 3), C5(1, 3), C6(1, 3);
	rowvec FreeVertex(1, 3);

	mat RepmatFreeVertex(9, 3);
	mat pho_Plus(EdgesTotal, 3), pho_Minus(EdgesTotal, 3);

	cube center_(9, 3, TrianglesTotal);
	cube pho_Plus_(9, 3, EdgesTotal), pho_Minus_(9, 3, EdgesTotal);

	for (int im = 0; im < TrianglesTotal; im++)
	{
		//组成一个 element 的三个点
		n1 = t_data(im, 0) - 1;  
		n2 = t_data(im, 1) - 1;
		n3 = t_data(im, 2) - 1;

		//一个 element 的三个点的坐标
		r1 = p_data.row(n1);     
		r2 = p_data.row(n2);
		r3 = p_data.row(n3);

		r12 = r2 - r1;
		r23 = r3 - r2;
		r13 = r3 - r1;

		M = Center.row(im);
		C1 = r1 + (1.0 / 3.0)*r12;
		C2 = r1 + (2.0 / 3.0)*r12;
		C3 = r2 + (1.0 / 3.0)*r23;
		C4 = r2 + (2.0 / 3.0)*r23;
		C5 = r1 + (1.0 / 3.0)*r13;
		C6 = r1 + (2.0 / 3.0)*r13;

		mat temp1(9, 3);
		temp1.row(0) = (1.0 / 3.0)*(C1 + C5 + r1);
		temp1.row(1) = (1.0 / 3.0)*(C1 + C2 + M);
		temp1.row(2) = (1.0 / 3.0)*(C2 + C3 + r2);
		temp1.row(3) = (1.0 / 3.0)*(C2 + C3 + M);
		temp1.row(4) = (1.0 / 3.0)*(C3 + C4 + M);
		temp1.row(5) = (1.0 / 3.0)*(C1 + C5 + M);
		temp1.row(6) = (1.0 / 3.0)*(C5 + C6 + M);
		temp1.row(7) = (1.0 / 3.0)*(C4 + C6 + M);
		temp1.row(8) = (1.0 / 3.0)*(C4 + C6 + r3);
		center_.slice(im) = temp1;
	}

	for (int im = 0; im < EdgesTotal; im++)
	{
		// 被定义为正的element     的编号
		NoPlus = TrianglePlus(im);
		//组成一个 element 的三个点
		n1 = t_data(NoPlus, 0) - 1;
		n2 = t_data(NoPlus, 1) - 1;
		n3 = t_data(NoPlus, 2) - 1;

		if (n1 != EdgeElement(im, 0) && n1 != EdgeElement(im, 1))
			NODE = n1;
		else if (n2 != EdgeElement(im, 0) && n2 != EdgeElement(im, 1))
			NODE = n2;
		else
			NODE = n3;

		FreeVertex = p_data.row(NODE);  //三角形公共边对应的另外一个点坐标
		for (int i = 0; i < 9; i++)
		{
			//mat(9,3)
			RepmatFreeVertex.row(i) = FreeVertex;
		}

		//非公共边点到重心的向量坐标 
		pho_Plus.row(im) = +Center.row(NoPlus) - FreeVertex;
		//九点划分后的 非公共边点到  小三角形重心的向量坐标
		pho_Plus_.slice(im) = +center_.slice(NoPlus) - RepmatFreeVertex;
	}

	for (int im = 0; im < EdgesTotal; im++)
	{
		NoMinus = TriangleMinus(im);
		n1 = t_data(NoMinus, 0) - 1;
		n2 = t_data(NoMinus, 1) - 1;
		n3 = t_data(NoMinus, 2) - 1;

		if (n1 != EdgeElement(im, 0) && n1 != EdgeElement(im, 1))
			NODE = n1;
		else if (n2 != EdgeElement(im, 0) && n2 != EdgeElement(im, 1))
			NODE = n2;
		else
			NODE = n3;

		FreeVertex = p_data.row(NODE);
		for (int i = 0; i < 9; i++)
		{
			RepmatFreeVertex.row(i) = FreeVertex;
		}

		pho_Minus.row(im) = -Center.row(NoMinus) + FreeVertex;
		pho_Minus_.slice(im) = -center_.slice(NoMinus) + RepmatFreeVertex;
	}

	Center_ = center_;
	Pho_Plus = pho_Plus;
	Pho_Minus = pho_Minus;
	Pho_Plus_ = pho_Plus_;
	Pho_Minus_ = pho_Minus_;
}
