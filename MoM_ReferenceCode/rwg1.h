#pragma once
/***********************************************************/
/*         �ú�����Ҫ������������ʷ����ݽ��д���          */
/*           �������εıߡ����������ε�Ԫ���б�            */
/*             �ţ�������������������ߵĳ���              */
/***********************************************************/
/*                                                         */
/*                   Author��Yanlin Xu                     */
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

int TrianglesTotal, EdgesTotal;        //�ֱ��ʾ������Ŀ���������Ƭ���Լ�RWG��������
mat t_data;                                       //�þ���洢Triangle�ڵ�����Ϣ
mat	p_data;										 //�þ���洢Triangle�ڵ�������Ϣ
mat EdgeElement;                           //�þ���洢EdgeElement�ڵ�����Ϣ
vec TrianglePlus;                            //�þ���洢EdgeElement��Ӧ���������α����Ϣ
vec TriangleMinus;                         //�þ���洢EdgeElement��Ӧ�ĸ������α����Ϣ     
vec EdgeLength;                              //�þ���洢EdgeElement�ĳ���
mat Center;                                      //�þ���洢Triangle��Ԫ����������
vec Area;                                          //������洢Triangle��Ԫ�����
extern string  t_filename, p_filename;
extern rowvec  FeedPoint;          //�����ѹԴ�������(only for r_s = 2)
extern rowvec  FeedPlus;            //��������������(only for r_s = 2)
extern rowvec  FeedEdge;           //��������ߵ�EdgeElement���
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
/*              Function��void edgenumber()                */
/***********************************************************/
/* 1.ͳ��EdgesTotal                                        */
/***********************************************************/
/*Ѱ�ұ�*/
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
/*             Function��void edgeelement()                */
/***********************************************************/
/* 1.ͳ��EdgesTotal                                        */
/* 2.��EdgeElement��ţ��������˽ڵ�                     */
/* 3.��TrianglePlus��TriangleMinus���                     */
/***********************************************************/
void edgeelement()
{
	int counter = 0;      //�ڴ����ڼ���EdgeElement����
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
	/*--------��һ�γ�������������������е�T junction---------*/
	int T_count = 0;
	mat edge__(EdgesTotal, 2);
	edge__.col(0) = edgeElement.col(1);
	edge__.col(1) = edgeElement.col(0);
	int common_edge_count = 0;
	rowvec common_edge_index(10000);
	rowvec remove(10000);
	for (int mm = 0; mm < EdgesTotal; mm++) {
		mat Edge_m = edgeElement.row(mm);      //edgeelement ���֮ǰ�ڴ洢ʱ�Ѿ���1
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
					rowvec trim = t_data.row(triangleMinus(common_edge_index(pp))) - 1;  ///���������⣡������
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
					rowvec trim = t_data.row(triangleMinus(common_edge_index(pp))) - 1;  ///���������⣡������
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
					rowvec trim = t_data.row(triangleMinus(common_edge_index(pp))) - 1;  ///���������⣡������
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
					rowvec trim = t_data.row(triangleMinus(common_edge_index(pp))) - 1;  ///���������⣡������
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
	cout << "\n" << "����⵽" << T_count << "��T junction," << "RWG����������" << EdgesTotal << "��������" << EdgesTotal - T_count << "��!" << endl;
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
/*              Function��void feededge()                  */
/***********************************************************/
/* 1.Ѱ�������                                            */
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
	cout << "\n" << "����⵽" << count << "������ߣ�" << endl;
}
/***********************************************************/
/*                 Function��float cross()                 */
/***********************************************************/
/* 1.�����������Ĳ��                                      */
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
/*              Function��float edgelength()               */
/***********************************************************/
/* 1.����EdgeElement�ĳ���                                 */
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
/*                 Function��float area()                  */
/***********************************************************/
/* 1.����Triangles�����                                   */
/* 2.����Triangles������                                   */
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
