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
extern vec EdgeLength;  //�����߳���
extern vec Vertex_pos;  // �������ε� �ǹ����ߵ� ���
extern vec Vertex_naga;

extern cube Center_;   //�ŵ���ֻ��ֺ��С����������

extern mat Pho_Pos;       //������element������ �ǹ����ߵ㵽���ĵ��������� 
extern mat Pho_Naga;      //���ڸ�element������ �ǹ����ߵ㵽���ĵ���������
extern cube Pho_Pos_;     //������element������ �ŵ㻮�ֺ� �ǹ����ߵ㵽 С���������ĵ���������
extern cube Pho_Naga_;    //���ڸ�element������ �ŵ㻮�ֺ� �ǹ����ߵ㵽 С���������ĵ���������

extern cx_mat Z;
extern rowvec Pol, dz;        //�������䲨�ļ�����������䷽��(only for r_s = 1)
extern double FeedVoltage;    //���������ѹ
extern rowvec FeedEdge;       //��������ߵ�EdgeElement���

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
	printf("\n ���������󷽳� ... \n");
	clock_t start = clock();
	Ii = solve(Z, Vv);          //ע����ĵĽ��� RWG���� ��ϵ��
	clock_t end = clock();
	cout << endl << "�ⷽ�̵�ʱ��Ϊ��" << (end - start) << "ms" << endl;
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
	//printf("\n ���������󷽳� ... \n");
	//clock_t start = clock();
	//Ii = solve(Z, Vv);          //ע����ĵĽ��� RWG���� ��ϵ��
	//clock_t end = clock();
	//cout << endl << "�ⷽ�̵�ʱ��Ϊ��" << (end - start) << "ms" << endl;
	V = Vv;
	//I = Ii;
}