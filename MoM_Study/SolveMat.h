#pragma once
#include<iostream>
#include<stdio.h>
#include<armadillo>
#include<math.h>
#include<complex>
#include<time.h>
#include<cmath>
#include"FillZmat.h"

using namespace arma;
using namespace std;

extern double k;
extern mat basis_list;
extern cube element_list;
extern int elem_num;
extern int rwg_num;
extern mat Center;
extern vec Area;

extern mat CommonEdgeP;
extern mat CommonEdgeN;
extern vec EdgeLength;  //公共边长度
extern vec Vertex_pos;  // 正三角形的 非公共边点 编号
extern vec Vertex_naga;

extern cube Center_;   //九点积分划分后的小三角形重心

extern mat Pho_Pos;       //对于正element三角形 非公共边点到重心的向量坐标 
extern mat Pho_Naga;      //对于负element三角形 非公共边点到重心的向量坐标
extern cube Pho_Pos_;     //对于正element三角形 九点划分后 非公共边点到 小三角形重心的向量坐标
extern cube Pho_Naga_;    //对于负element三角形 九点划分后 非公共边点到 小三角形重心的向量坐标

extern cx_mat Z;
extern rowvec Pol, dz;        //定义入射波的极化方向和入射方向(only for r_s = 1)
extern double FeedVoltage;    //定义馈电电压
extern rowvec FeedEdge;       //定义馈电边的EdgeElement编号

rowvec kv(1, 3);
cx_colvec I, V;
double P_rad = 0;

void SolveMat() {
	kv = k * dz;
	FillmatZ();
	double ScalarTemp = 0.0;
	cx_rowvec Einc_Pos(1, 3), Einc_Nega(1, 3); 
	cx_colvec Vv(rwg_num, 1), Ii(rwg_num, 1);
	complex<double>EPho_Pos(0, 0), EPho_Nega(0, 0);

	for (int i = 0; i < rwg_num; i++) {
		ScalarTemp = as_scalar(kv * strans(Center.row(basis_list(i, 0))));
		Einc_Pos = Pol * exp(complex<double>(0, -ScalarTemp));

		ScalarTemp = as_scalar(kv * strans(Center.row(basis_list(i, 1))));
		Einc_Nega = Pol * exp(complex<double>(0, -ScalarTemp));

		EPho_Pos = as_scalar(Einc_Pos * strans(Pho_Pos.row(i)));
		EPho_Nega = as_scalar(Einc_Nega * strans(Pho_Naga.row(i)));
		Vv = EdgeLength(i) * (EPho_Pos / 2.0 + EPho_Nega / 2.0);
	}
	printf("\n 正在求解矩阵方程 ... \n");
	clock_t start = clock();
	Ii = solve(Z, Vv);          //注意求的的解是 RWG函数 的系数
	clock_t end = clock();
	cout << endl << "解方程的时间为：" << (end - start) << "ms" << endl;
	V = Vv;
	I = Ii;
}

void SolveMat_V() {
	kv = k * dz;
	//FillmatZ();
	double ScalarTemp = 0.0;
	cx_rowvec Einc_Pos(1, 3), Einc_Nega(1, 3);
	cx_colvec Vv(rwg_num, 1), Ii(rwg_num, 1);
	complex<double>EPho_Pos(0, 0), EPho_Nega(0, 0);

	for (int i = 0; i < rwg_num; i++) {
		ScalarTemp = as_scalar(kv * strans(Center.row(basis_list(i, 0))));
		Einc_Pos = Pol * exp(complex<double>(0, -ScalarTemp));

		ScalarTemp = as_scalar(kv * strans(Center.row(basis_list(i, 1))));
		Einc_Nega = Pol * exp(complex<double>(0, -ScalarTemp));

		EPho_Pos = as_scalar(Einc_Pos * strans(Pho_Pos.row(i)));
		EPho_Nega = as_scalar(Einc_Nega * strans(Pho_Naga.row(i)));
		Vv = EdgeLength(i) * (EPho_Pos / 2.0 + EPho_Nega / 2.0);
	}
	//printf("\n 正在求解矩阵方程 ... \n");
	//clock_t start = clock();
	//Ii = solve(Z, Vv);          //注意求的的解是 RWG函数 的系数
	//clock_t end = clock();
	//cout << endl << "解方程的时间为：" << (end - start) << "ms" << endl;
	V = Vv;
	//I = Ii;
}