#pragma once
/***********************************************************/
/*         该函数主要用来对输入的剖分数据进行处理          */
/*           对三角形的边、正负三角形单元进行编            */
/*             号，并计算三角形面积及边的长度              */
/***********************************************************/
/*                                                         */
/*                   Author：Yanlin Xu                     */
/*                                                         */
/***********************************************************/
#include "rwg1.h"
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "armadillo"

using namespace arma;
using namespace std;

int TrianglesTotal, EdgesTotal;        //分别表示待分析目标的三角面片数以及RWG函数数量
mat t_data;                                       //该矩阵存储Triangle节点编号信息
mat	p_data;										 //该矩阵存储Triangle节点坐标信息
mat EdgeElement;                           //该矩阵存储EdgeElement节点编号信息
vec TrianglePlus;                            //该矩阵存储EdgeElement对应的正三角形编号信息
vec TriangleMinus;                         //该矩阵存储EdgeElement对应的负三角形编号信息     
vec EdgeLength;                              //该矩阵存储EdgeElement的长度
mat Center;                                      //该矩阵存储Triangle单元的重心坐标
vec Area;                                          //该数组存储Triangle单元的面积
extern string  t_filename, p_filename;
extern rowvec  FeedPoint;          //定义电压源的馈电点(only for r_s = 2)
extern rowvec  FeedPlus;            //定义馈电点的正极(only for r_s = 2)
extern rowvec  FeedEdge;           //定义馈电边的EdgeElement编号
extern mat FeedMat;
extern mat FeedPlusMat;

void rwg1_main()
{
	mat Tdata, Pdata;
	TrianglesTotal = 0;
	EdgesTotal = 0;
	Tdata.load(t_filename, raw_ascii);
	Pdata.load(p_filename, raw_ascii);
	t_data = Tdata;
	p_data = Pdata;
	TrianglesTotal = t_data.n_rows;
	area();
	edgenumber();
	edgeelement();
	edgelength();
	cout << "\n" << "TrianglesTotal: " << TrianglesTotal << endl;
	cout << "\n" << "EdgesTotal: " << EdgesTotal << endl;
	feededge();
	//EdgeLength.save("EdgeLength.txt", raw_ascii);
	//EdgeElement.save("EdgeElement.txt",raw_ascii);
}
/***********************************************************/
/*              Function：void edgenumber()                */
/***********************************************************/
/* 1.统计EdgesTotal                                        */
/***********************************************************/
/*寻找边*/
void edgenumber()
{
	int N[3] = { 0 }, M[3] = { 0 };
	for (int ia = 0; ia < TrianglesTotal; ia++)
	{
		for (int i = 0; i < 3; i++)
			N[i] = t_data(ia, i) - 1;
		for (int ib = ia + 1; ib < TrianglesTotal; ib++)
		{
			for (int j = 0; j < 3; j++)
				M[j] = t_data(ib, j) - 1;
			int samenode = 0;
			for (int ic = 0; ic < 3; ic++)
			{
				for (int id = 0; id < 3; id++)
				{
					if (N[ic] == M[id])
						samenode = samenode + 1;
				}
			}
			if (samenode == 2)
			{
				EdgesTotal++;
			}
		}
	}
}
/***********************************************************/
/*             Function：void edgeelement()                */
/***********************************************************/
/* 1.统计EdgesTotal                                        */
/* 2.对EdgeElement编号，并绑定两端节点                     */
/* 3.对TrianglePlus、TriangleMinus编号                     */
/***********************************************************/
void edgeelement()
{
	int counter = 0;      //在此用于记数EdgeElement数量
	mat edgeElement(EdgesTotal, 2);
	vec trianglePlus(EdgesTotal);
	vec triangleMinus(EdgesTotal);
	int N[3] = { 0 }, M[3] = { 0 };
	for (int ia = 0; ia < TrianglesTotal; ia++)
	{
		for (int i = 0; i < 3; i++)
			N[i] = t_data(ia, i) - 1;
		for (int ib = ia + 1; ib < TrianglesTotal; ib++)
		{
			for (int j = 0; j < 3; j++)
				M[j] = t_data(ib, j) - 1;
			int samenode = 0;
			for (int ic = 0; ic < 3; ic++)
			{
				for (int id = 0; id < 3; id++)
				{
					if (N[ic] == M[id])
						samenode = samenode + 1;
				}
			}
			if (samenode == 2)
			{
				counter++;
				samenode = 2;
				for (int ic = 0; ic < 3; ic++)
				{
					for (int id = 0; id < 3; id++)
					{
						if (N[ic] == M[id])
						{
							samenode = samenode - 1;
							edgeElement(counter - 1, samenode) = M[id];
						}
					}
				}
				trianglePlus(counter - 1) = ia;
				triangleMinus(counter - 1) = ib;
			}
		}
	}
	/*--------这一段程序用来检测三角网格中的T junction---------*/
	int T_count = 0;
	mat edge__(EdgesTotal, 2);
	edge__.col(0) = edgeElement.col(1);
	edge__.col(1) = edgeElement.col(0);
	int common_edge_count = 0;
	rowvec common_edge_index(10000);
	rowvec remove(10000);
	for (int mm = 0; mm < EdgesTotal; mm++) {
		mat Edge_m = edgeElement.row(mm);      //edgeelement 编号之前在存储时已经减1
		mat Edge_m_center = (p_data.row(Edge_m(0)) + p_data.row(Edge_m(1))) / 2;
		common_edge_count = 1;
		common_edge_index.zeros();
		common_edge_index(common_edge_count - 1) = mm;
		for (int nn = mm + 1; nn < EdgesTotal; nn++) {
			mat ind1 = edgeElement.row(nn) - Edge_m;
			mat ind2 = edge__.row(nn) - Edge_m;
			int Ind1 = abs(ind1(0)) + abs(ind1(1));
			int Ind2 = abs(ind2(0)) + abs(ind2(1));
			if (Ind1*Ind2 == 0) {
				common_edge_count++;
				common_edge_index(common_edge_count - 1) = nn;
			}
		}

		if (common_edge_count == 3) {
			T_count++;
			/*if(norm(Edge_m_center-FeedPoint,2) == 0){
				for(int pp = 0;pp<3;pp++){
					rowvec trip = t_data.row(trianglePlus(common_edge_index(pp))) - 1;
					rowvec trim = t_data.row(triangleMinus(common_edge_index(pp))) - 1;  ///这里有问题！！！！
					int vertexp = accu(trip) - accu(Edge_m);
					int vertexm = accu(trim) - accu(Edge_m);
					if( norm(p_data.row(vertexp)-FeedPlus,2) == 0 || norm(p_data.row(vertexm)-FeedPlus,2) == 0 )
						continue;
					else{
						remove(T_count-1) = common_edge_index(pp);
						break;
					}
				}
			}*/
			if (norm(Edge_m_center - FeedMat.row(0), 2) == 0) {
				FeedPlus = FeedPlusMat.row(0);
				for (int pp = 0; pp < 3; pp++) {
					rowvec trip = t_data.row(trianglePlus(common_edge_index(pp))) - 1;
					rowvec trim = t_data.row(triangleMinus(common_edge_index(pp))) - 1;  ///这里有问题！！！！
					int vertexp = accu(trip) - accu(Edge_m);
					int vertexm = accu(trim) - accu(Edge_m);
					if (norm(p_data.row(vertexp) - FeedPlus, 2) == 0 || norm(p_data.row(vertexm) - FeedPlus, 2) == 0)
						continue;
					else {
						remove(T_count - 1) = common_edge_index(pp);
						break;
					}
				}
			}
			else if (norm(Edge_m_center - FeedMat.row(1), 2) == 0) {
				FeedPlus = FeedPlusMat.row(1);
				for (int pp = 0; pp < 3; pp++) {
					rowvec trip = t_data.row(trianglePlus(common_edge_index(pp))) - 1;
					rowvec trim = t_data.row(triangleMinus(common_edge_index(pp))) - 1;  ///这里有问题！！！！
					int vertexp = accu(trip) - accu(Edge_m);
					int vertexm = accu(trim) - accu(Edge_m);
					if (norm(p_data.row(vertexp) - FeedPlus, 2) == 0 || norm(p_data.row(vertexm) - FeedPlus, 2) == 0)
						continue;
					else {
						remove(T_count - 1) = common_edge_index(pp);
						break;
					}
				}
			}
			else if (norm(Edge_m_center - FeedMat.row(2), 2) == 0) {
				FeedPlus = FeedPlusMat.row(2);
				for (int pp = 0; pp < 3; pp++) {
					rowvec trip = t_data.row(trianglePlus(common_edge_index(pp))) - 1;
					rowvec trim = t_data.row(triangleMinus(common_edge_index(pp))) - 1;  ///这里有问题！！！！
					int vertexp = accu(trip) - accu(Edge_m);
					int vertexm = accu(trim) - accu(Edge_m);
					if (norm(p_data.row(vertexp) - FeedPlus, 2) == 0 || norm(p_data.row(vertexm) - FeedPlus, 2) == 0)
						continue;
					else {
						remove(T_count - 1) = common_edge_index(pp);
						break;
					}
				}
			}
			else
				remove(T_count - 1) = common_edge_index(0);
		}
	}
	cout << "\n" << "共检测到" << T_count << "个T junction," << "RWG函数数量从" << EdgesTotal << "个缩减到" << EdgesTotal - T_count << "个!" << endl;
	if (T_count > 0) {
		rowvec Remove = sort(remove.cols(0, T_count - 1), 1);  //sort: from larger to smaller 
		for (int mm = 0; mm < T_count; mm++) {
			edgeElement.shed_row(Remove(mm));
			trianglePlus.shed_row(Remove(mm));
			triangleMinus.shed_row(Remove(mm));
		}
	}
	EdgesTotal = EdgesTotal - T_count;
	EdgeElement = edgeElement;
	TrianglePlus = trianglePlus;
	TriangleMinus = triangleMinus;
}
/***********************************************************/
/*              Function：void feededge()                  */
/***********************************************************/
/* 1.寻找馈电边                                            */
/***********************************************************/
void feededge()
{
	rowvec fe(10000);
	int count = 0;
	for (int mm = 0; mm < EdgesTotal; mm++) {
		mat Edge_m = EdgeElement.row(mm);
		mat Edge_m_center = (p_data.row(Edge_m(0)) + p_data.row(Edge_m(1))) / 2;
		for (int nn = 0; nn < FeedMat.n_rows; nn++) {
			if (norm(Edge_m_center - FeedMat.row(nn), 2) == 0) {
				fe(count) = mm;
				count++;
				break;
			}
		}
	}
	if (count > 0)
		FeedEdge = fe.cols(0, count - 1);
	cout << "\n" << "共检测到" << count << "个馈电边！" << endl;
}
/***********************************************************/
/*                 Function：float cross()                 */
/***********************************************************/
/* 1.求两个向量的差乘                                      */
/***********************************************************/
float * cross(float coor1[3], float coor2[3])
{
	float * coor3 = new float[3];
	coor3[0] = coor1[1] * coor2[2] - coor1[2] * coor2[1];
	coor3[1] = coor1[2] * coor2[0] - coor1[0] * coor2[2];
	coor3[2] = coor1[0] * coor2[1] - coor1[1] * coor2[0];
	return coor3;
}
/***********************************************************/
/*              Function：float edgelength()               */
/***********************************************************/
/* 1.计算EdgeElement的长度                                 */
/***********************************************************/
void edgelength()
{
	vec edgeLength(EdgesTotal);
	rowvec vector1(1, 3), vector2(1, 3), vector3(1, 3);
	for (int ia = 0; ia < EdgesTotal; ia++)
	{

		vector1 = p_data.row(EdgeElement(ia, 0));
		vector2 = p_data.row(EdgeElement(ia, 1));
		vector3 = vector1 - vector2;
		edgeLength(ia) = sqrt(as_scalar(vector3*strans(vector3)));
	}
	EdgeLength = edgeLength;
}
/***********************************************************/
/*                 Function：float area()                  */
/***********************************************************/
/* 1.计算Triangles的面积                                   */
/* 2.计算Triangles的重心                                   */
/***********************************************************/
void area()
{
	vec area(TrianglesTotal);
	mat center(TrianglesTotal, 3);
	float test = 0, test1 = 0;
	float vector1[3] = { 0 }, vector2[3] = { 0 };
	float * vector3 = new float[3];
	int N[3] = { 0 };
	for (int ia = 0; ia < TrianglesTotal; ia++)
	{
		for (int i = 0; i < 3; i++)
			N[i] = t_data(ia, i) - 1;
		for (int j = 0; j < 3; j++)
		{
			vector1[j] = p_data(N[0], j) - p_data(N[1], j);
			vector2[j] = p_data(N[2], j) - p_data(N[1], j);
		}
		vector3 = cross(vector1, vector2);
		area(ia) = 0.5*sqrt(vector3[0] * vector3[0] + vector3[1] * vector3[1] + vector3[2] * vector3[2]);
		center(ia, 0) = (1.0 / 3.0)*(p_data(N[0], 0) + p_data(N[1], 0) + p_data(N[2], 0));
		center(ia, 1) = (1.0 / 3.0)*(p_data(N[0], 1) + p_data(N[1], 1) + p_data(N[2], 1));
		center(ia, 2) = (1.0 / 3.0)*(p_data(N[0], 2) + p_data(N[1], 2) + p_data(N[2], 2));
	}
	Area = area;
	Center = center;
}
