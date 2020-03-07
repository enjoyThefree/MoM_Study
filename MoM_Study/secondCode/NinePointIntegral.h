#pragma once
#include<iostream>
#include<string>
#include<armadillo>
#include<math.h>
#include"InputData.h"

using namespace std;
using namespace arma;

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

cube Center_;   //九点积分划分后的小三角形重心

mat Pho_Pos;       //对于正element三角形 非公共边点到重心的向量坐标 
mat Pho_Naga;      //对于负element三角形 非公共边点到重心的向量坐标
cube Pho_Pos_;     //对于正element三角形 九点划分后 非公共边点到 小三角形重心的向量坐标
cube Pho_Naga_;    //对于负element三角形 九点划分后 非公共边点到 小三角形重心的向量坐标

void NinepointIntegral() {

	InputData();
	cube center_(9, 3, elem_num);
	rowvec r1(1, 3), r2(1, 3), r3(1, 3);
	rowvec r12(1, 3), r23(1, 3), r13(1, 3);
	rowvec M(1, 3), c1(1, 3), c2(1, 3), c3(1, 3), c4(1, 3), c5(1, 3), c6(1, 3);

	int NoVertexP = 0; int NoVertexN = 0;   //非公共边点 编号
	
	rowvec VertexP(1, 3);   // 正三角形的 非公共边点 坐标信息
	rowvec VertexN(1, 3);   // 负三角形的 非公共边点 坐标信息
	mat FreeVertex_Pmat(9, 3);    //九点划分后,正三角形的 扩展的 非公共边点 坐标信息
	mat FreeVertex_Nmat(9, 3);    //九点划分后,正三角形的 扩展的 非公共边点 坐标信息
	mat pho_pos(rwg_num, 3);     
	mat pho_nega(rwg_num, 3);    

	cube pho_pos_(9, 3, rwg_num);   
	cube pho_nega_(9, 3, rwg_num);

	for (int i = 0; i < elem_num; i++) {
		//组成一个 element 的三个点的坐标；
		r1 = element_list.slice(i).row(0);
		r2 = element_list.slice(i).row(1);
		r3 = element_list.slice(i).row(2);

		//三角形三条边的向量
		r12 = r2 - r1;
		r23 = r3 - r2;
		r13 = r3 - r1;

		// element 的重心
		M = Center.row(i);

		// element 三角形边上的三等分点
		c1 = r1 + (1 / 3.0) * r12;
		c2 = r1 + (2 / 3.0) * r12;

		c3 = r2 + (1 / 3.0) * r23;
		c4 = r2 + (2 / 3.0) * r23;

		c5 = r1 + (1 / 3.0) * r13;
		c5 = r1 + (2 / 3.0) * r13;

		//计算每个小三角形的重心
		mat temp(9,3);
		temp.row(0) = (1 / 3.0) * (c1 + c3 + r1);
		temp.row(1) = (1 / 3.0) * (c1 + c2 + M);
		temp.row(2) = (1 / 3.0) * (c2 + c3 + r2);
		temp.row(3) = (1 / 3.0) * (c2 + c3 + M);
		temp.row(4) = (1 / 3.0) * (c3 + c4 + M);
		temp.row(5) = (1 / 3.0) * (c1 + c5 + M);
		temp.row(6) = (1 / 3.0) * (c5 + c6 + M);
		temp.row(7) = (1 / 3.0) * (c4 + c6 + M);
		temp.row(8) = (1 / 3.0) * (c4 + c6 + r3);

		center_.slice(i) = temp;
	}

	for (int i = 0; i < rwg_num; i++) {
		// 取出 被定义为正、负rwg的element 的编号
		int NoPos = basis_list(i, 0);
		int NoNega = basis_list(i, 1);
		// 取出 正、负三角形的 非公共边点 编号
		NoVertexP = Vertex_pos(i);
		NoVertexN = Vertex_naga(i);

		VertexP = element_list.slice(NoPos).row(NoVertexP);
		VertexN = element_list.slice(NoNega).row(NoVertexN);

		for (int i = 0; i < 9; i++) {
			FreeVertex_Pmat.row(i) = VertexP;
			FreeVertex_Nmat.row(i) = VertexN;
		}

		//非公共边点到重心的向量坐标 
		pho_pos.row(i) = Center.row(NoPos) - VertexP;
		pho_pos_.slice(i) = center_.slice(NoPos) - FreeVertex_Pmat;

		pho_nega.row(i) = Center.row(NoNega) - VertexN;
		pho_nega_.slice(i) = center_.slice(NoNega) - FreeVertex_Nmat;

	}

	Center_ = center_;
	Pho_Pos = pho_pos;
	Pho_Pos_ = pho_pos_;
	Pho_Naga = pho_nega;
	Pho_Naga_ = pho_nega_;
}