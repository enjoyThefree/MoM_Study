#pragma once
/***********************************************************/
/*                该函数是填充阻抗矩阵元素                 */
/*                     定义入射波参数                      */
/***********************************************************/
/*                                                         */
/*                   Author：Yanlin Xu                     */
/*                                                         */
/***********************************************************/
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <complex>
#include <math.h>
#include <time.h>
#include <windows.h>
#include <cmath>
#include "rwg3.h"
#include "armadillo"

using namespace arma;
using namespace std;

extern int      TrianglesTotal, EdgesTotal;
extern mat      t_data;
extern mat      p_data;
extern mat      EdgeElement;
extern vec      TrianglePlus;     //该矩阵存储EdgeElement对应的正三角形编号信息
extern vec      TriangleMinus;    //该矩阵存储EdgeElement对应的负三角形编号信息
extern vec      EdgeLength;
extern mat      Center;
extern vec      Area;
extern mat      Pho_Plus;       //  正三角形 非公共边对应点 到 重心 的矢量 扩展 9*3*rwg_num
extern mat      Pho_Minus;      //  负三角形 非公共边对应点 到 重心 的矢量 扩展 9*3*rwg_num
extern cube     Center_;
extern cube     Pho_Plus_;
extern cube     Pho_Minus_;
extern float    f;
extern float    c;
extern float    lamda;
extern float    eta;
extern float    omega;
extern double   k;    //波数
extern int      SimpleGF;
extern double   delta;
extern cx_rowvec g_base;
extern mat      LoadMat;

extern double mu = 4 * pi * 1e-7;
extern double epsilon = 8.854187817 * 1e-12;

float constant1 = mu / (4 * pi);
float factor = 1.0 / 9.0;
float Z_fill_time = 0.0;
#define pi2     6.283184

cx_mat Z;

void rwg3_main()
{


	int GFindex = 0;
	double Index_temp = 0;


	complex<float>  K(0, k);
	complex<float>  constant2(0, -1 / (4 * pi*omega*epsilon));

	rwg2_main();
	cx_colvec Factor_A(EdgesTotal, 1), Factor_Fi(EdgesTotal, 1);
	cx_colvec Fi(EdgesTotal, 1), ZF(EdgesTotal, 1), Z1(EdgesTotal, 1);
	cx_colvec A(EdgesTotal, 1);

	rowvec Plus, Minus;    // 正负三角形编号

	cube Pho_P(9, 3, EdgesTotal), Pho_M(9, 3, EdgesTotal);    //   非公共边对应点 到 重心 的矢量 扩展 9*3*rwg_num
	cube RepmatCenter(9, 3, TrianglesTotal);     // element 重心 扩展
	cube RepmatCenter_(9, 3, TrianglesTotal);    // 九点划分后的重心
	cube RP(9, 3, EdgesTotal);                   //对于正element三角形 九点划分后 非公共边点到 小三角形重心的向量坐标

	mat D(9, 3);   // element原三角形重心 减去 九点划分后的小三角形重心 的向量
	mat R(TrianglesTotal, 9);    // R = |r - r'|
	cx_mat g(TrianglesTotal, 9);
	cx_mat gP(EdgesTotal, 9), gM(EdgesTotal, 9);
	cx_mat z(EdgesTotal, EdgesTotal);

	for (int im = 0; im < EdgesTotal; im++)
	{
		Factor_A(im) = complex<float>(0, factor*(omega*EdgeLength(im) / 4.0)*constant1);
		Factor_Fi(im) = complex<float>(factor*EdgeLength(im), 0)*constant2;

		mat temp_P(9, 3), temp_M(9, 3);
		for (int i = 0; i < 9; i++)
		{
			temp_P.row(i) = Pho_Plus.row(im);
			temp_M.row(i) = Pho_Minus.row(im);
		}
		Pho_P.slice(im) = temp_P;   //第m个element的   非公共边对应点 到 重心 的矢量 扩展 9*3*rwg_num
		Pho_M.slice(im) = temp_M;
	}

	cout << "\n" << "正在计算阻抗矩阵，请等待..." << endl;
	clock_t start = clock();
	//double deleteTime = 0;
	z = zeros<cx_mat>(EdgesTotal, EdgesTotal);
	for (int ip = 0; ip < TrianglesTotal; ip++)   // ip element编号
	{
		//clock_t dT_start = clock();
		int length_Plus = 0, length_Minus = 0;
		// 在RWG中 正负三角形element编号 所在位置
		for (int im = 0; im < EdgesTotal; im++)  //rwg
		{
			if (TrianglePlus(im) == ip)
				length_Plus++;
			if (TriangleMinus(im) == ip)
				length_Minus++;
		}
		if (length_Plus > 0)
		{
			int io = 0;
			rowvec plus(length_Plus);
			for (int im = 0; im < EdgesTotal; im++)
			{
				if (TrianglePlus(im) == ip)
				{
					plus(io) = im;
					io++;
				}
			}
			Plus = plus;
		}
		if (length_Minus > 0)
		{
			int io = 0;
			rowvec minus(length_Minus);
			for (int im = 0; im < EdgesTotal; im++)
			{
				if (TriangleMinus(im) == ip)
				{
					minus(io) = im;
					io++;
				}
			}
			Minus = minus;
		}

		for (int im = 0; im < TrianglesTotal; im++)       // im element 序号
		{
			RepmatCenter_.slice(im) = Center_.slice(ip);

			mat temp(9, 3);
			for (int in = 0; in < 9; in++)
			{
				temp.row(in) = Center.row(im);
			}
			RepmatCenter.slice(im) = temp;

			D = RepmatCenter.slice(im) - RepmatCenter_.slice(im);   //
			for (int i = 0; i < 9; i++)
			{
				R(im, i) = sqrt(as_scalar(D.row(i)*strans(D.row(i))));
				if (SimpleGF == 0) {
					//g(im,i) = exp(-K*complex<float>(R(im,i),0))/complex<float>(R(im,i),0);
					double kR_temp = k * R(im, i);
					g(im, i) = exp(complex<double>(0, -kR_temp)) / R(im, i);

				}
				else if (SimpleGF == 1) {
					double kR_temp = k * R(im, i);
					Index_temp = (kR_temp - pi2 * floor(kR_temp / pi2)) / delta;
					GFindex = floor(Index_temp + 0.5);
					g(im, i) = g_base(GFindex) / R(im, i);
				}
			}
		}

		for (int in = 0; in < EdgesTotal; in++)
		{
			gP.row(in) = g.row(TrianglePlus(in));
			gM.row(in) = g.row(TriangleMinus(in));
		}
		//clock_t dT_end = clock();
		//deleteTime+=(dT_end-dT_start);

		Fi = sum(gP, 1) - sum(gM, 1);     //sum(A,dim),dim=0表示列相加，dim=1表示行相加
		ZF = Factor_Fi % Fi;

		int nn = 0;
		if (length_Plus > 0)
		{
			for (int im = 0; im < length_Plus; im++)
			{
				nn = Plus(im);
				for (int in = 0; in < EdgesTotal; in++)
				{
					RP.slice(in) = Pho_Plus_.slice(nn);
				}
				cube temp_p(9, 3, EdgesTotal), temp_m(9, 3, EdgesTotal);
				mat temp_P(EdgesTotal, 9), temp_M(EdgesTotal, 9);
				temp_p = RP % Pho_P;
				temp_m = RP % Pho_M;
				for (int i = 0; i < EdgesTotal; i++)
				{
					temp_P.row(i) = strans(sum(temp_p.slice(i), 1));    // 将一个element的 9个小三角形的 三个坐标相加；
					temp_M.row(i) = strans(sum(temp_m.slice(i), 1));
				}
				A = sum(gP % temp_P, 1) + sum(gM%temp_M, 1);
				Z1 = Factor_A % A;
				z.col(nn) = z.col(nn) + EdgeLength(nn)*(Z1 + ZF);
			}
		}

		nn = 0;
		if (length_Minus > 0)
		{
			for (int im = 0; im < length_Minus; im++)
			{
				nn = Minus(im);
				for (int in = 0; in < EdgesTotal; in++)
				{
					RP.slice(in) = Pho_Minus_.slice(nn);
				}
				cube temp_p(9, 3, EdgesTotal), temp_m(9, 3, EdgesTotal);
				mat temp_P(EdgesTotal, 9), temp_M(EdgesTotal, 9);
				temp_p = RP % Pho_P;
				temp_m = RP % Pho_M;
				for (int i = 0; i < EdgesTotal; i++)
				{
					temp_P.row(i) = strans(sum(temp_p.slice(i), 1));
					temp_M.row(i) = strans(sum(temp_m.slice(i), 1));
				}
				A = sum(gP%temp_P, 1) + sum(gM%temp_M, 1);
				Z1 = Factor_A % A;
				z.col(nn) = z.col(nn) + EdgeLength(nn)*(Z1 - ZF);
			}
		}
	}
	//-----寻找加载位置对应的RWG函数
	for (int im = 0; im < EdgesTotal; im++) {
		rowvec temp_node = EdgeElement.row(im);
		rowvec temp_center = (p_data.row(temp_node(0)) + p_data.row(temp_node(1))) / 2.0;
		if (norm(temp_center - LoadMat) < 1e-8)
		{
			z(im, im) = z(im, im) + EdgeLength(im)*EdgeLength(im) * 50;//complex<double>(0,omega*50.0);
			cout << "Find load" << endl;
		}

	}
	clock_t end = clock();
	//Z_fill_time = (end-start-deleteTime)/1000.0;
	cout << "\n" << "阻抗矩阵填充时间： " << (end - start) / 1000.0 << "s" << endl;
	//cout << deleteTime/1000 << "\t" << Z_fill_time << endl;
	Z = z;
}
