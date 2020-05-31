#pragma once
/*            ���ڼ�������ʵ���Ͼ���Ŀ����涨���RWG���������䲨���໥���ã�
                �ʶ��������Ŀ��������ʷ���Ϣ���֮�⣬�������䲨�ķ��򡢼�����ʽ��ء�
                ���⣬���������ı������Iʵ���ϲ�����Ŀ��ı�������ܶȣ���ֻ��RWG������ϵ����
              ���ߺ���һ������ϸ������ϵı�������ܶȡ�*/
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <complex>
#include <math.h>
#include <cmath>
#include "rwg4.h"
#include "armadillo"

using namespace arma;
using namespace std;

extern double k;
extern int    TrianglesTotal, EdgesTotal;
extern mat    t_data;         //�þ���洢Triangle�ڵ�����Ϣ
extern mat    p_data;         //�þ���洢Triangle�ڵ�������Ϣ
extern mat    EdgeElement;    //�þ���洢EdgeElement�ڵ�����Ϣ
extern vec    TrianglePlus;
extern vec    TriangleMinus;
extern vec    EdgeLength;
extern mat    Center;         //�þ���洢Triangle��Ԫ����������
extern vec    Area;           //������洢Triangle��Ԫ�����
extern mat    Pho_Plus;
extern mat    Pho_Minus;
extern cube   Center_;
extern cube   Pho_Plus_;
extern cube   Pho_Minus_;
extern cx_mat Z;
extern rowvec Pol, dz;
extern float  FeedVoltage;
extern rowvec FeedEdge;

rowvec kv(1, 3);
cx_colvec I, V;
float P_rad = 0;
complex<double> Zin;

void rwg4_main()
{
	kv = k * dz;
	rwg3_main();
	double ScalarProduct = 0.0;
	cx_rowvec EmPlus(1, 3), EmMinus(1, 3);
	cx_colvec Vv(EdgesTotal, 1), Ii(EdgesTotal);
	complex<double> ScalarPlus(0, 0), ScalarMinus(0, 0);
	for (int im = 0; im < EdgesTotal; im++)
	{
		ScalarProduct = as_scalar(kv*strans(Center.row(TrianglePlus(im))));
		EmPlus = Pol * exp(complex<double>(0, -ScalarProduct));
		ScalarProduct = as_scalar(kv*strans(Center.row(TriangleMinus(im))));
		EmMinus = Pol * exp(complex<double>(0, -ScalarProduct));
		ScalarPlus = as_scalar(EmPlus*strans(Pho_Plus.row(im)));
		ScalarMinus = as_scalar(EmMinus*strans(Pho_Minus.row(im)));
		Vv(im) = EdgeLength(im)*(ScalarPlus / 2.0 + ScalarMinus / 2.0);
	}
	cout << "\n" << "���������󷽳̣���ȴ�..." << endl;
	clock_t start = clock();
	Ii = solve(Z, Vv);
	clock_t end = clock();
	cout << "\n" << "���󷽳����ʱ�䣺 " << (end - start) / 1000.0 << "s" << endl;
	V = Vv;
	I = Ii;
}

void rwg4_monoV()
{
	kv = k * dz;
	double ScalarProduct = 0.0;
	cx_rowvec EmPlus(1, 3), EmMinus(1, 3);
	cx_colvec Vv(EdgesTotal, 1), Ii(EdgesTotal);
	complex<double> ScalarPlus(0, 0), ScalarMinus(0, 0);
	for (int im = 0; im < EdgesTotal; im++)
	{
		ScalarProduct = as_scalar(kv*strans(Center.row(TrianglePlus(im))));
		EmPlus = Pol * exp(complex<double>(0, -ScalarProduct));
		ScalarProduct = as_scalar(kv*strans(Center.row(TriangleMinus(im))));
		EmMinus = Pol * exp(complex<double>(0, -ScalarProduct));
		ScalarPlus = as_scalar(EmPlus*strans(Pho_Plus.row(im)));
		ScalarMinus = as_scalar(EmMinus*strans(Pho_Minus.row(im)));
		Vv(im) = EdgeLength(im)*(ScalarPlus / 2.0 + ScalarMinus / 2.0);
	}
	V = Vv;
}

void rwg4_radiate_V()
{
	rwg3_main();
	cx_colvec Vv(EdgesTotal, 1), Ii(EdgesTotal);
	Vv.zeros();
	for (int m = 0; m < FeedEdge.n_elem; m++)
		Vv(FeedEdge(m)) = FeedVoltage * EdgeLength(FeedEdge(m));
	cout << "\n" << "���������󷽳̣���ȴ�..." << endl;
	clock_t start = clock();
	Ii = solve(Z, Vv);
	clock_t end = clock();
	cout << "\n" << "���󷽳����ʱ�䣺 " << (end - start) / 1000.0 << "s" << endl;
	V = Vv;
	I = Ii;
	cx_rowvec VandI(1, 2);// Vgap = 0, Igap = 0;
	VandI.zeros();
	for (int m = 0; m < FeedEdge.n_elem; m++) {
		VandI(0) = VandI(0) + Vv(FeedEdge(m)) / EdgeLength(FeedEdge(m));
		VandI(1) = VandI(1) + Ii(FeedEdge(m))*EdgeLength(FeedEdge(m));
	}
	P_rad = 0.5*real(VandI(1)*conj(VandI(0))) / FeedEdge.n_elem;
	Zin = VandI(0) / VandI(1) / complex<double>(FeedEdge.n_elem, 0);
	//Zin = complex<double>(FeedVoltage,0)/VandI(1);
	cout << "\n" << "Zin: " << Zin << endl;
	//I.save("I_rv.txt",raw_ascii);
	//V.save("V_r.txt",raw_ascii);
	/*cx_rowvec VandI(1,6);// Vgap = 0, Igap = 0;
	VandI.zeros();
	for(int m=0;m<2;m++){
		VandI(0) = VandI(0) + Vv(FeedEdge(m))/EdgeLength(FeedEdge(m));
		VandI(1) = VandI(1) + Ii(FeedEdge(m))*EdgeLength(FeedEdge(m));
		VandI(2) = VandI(2) + Vv(FeedEdge(m+2))/EdgeLength(FeedEdge(m+2));
		VandI(3) = VandI(3) + Ii(FeedEdge(m+2))*EdgeLength(FeedEdge(m+2));
		VandI(4) = VandI(4) + Vv(FeedEdge(m+4))/EdgeLength(FeedEdge(m+4));
		VandI(5) = VandI(5) + Ii(FeedEdge(m+4))*EdgeLength(FeedEdge(m+4));
	}
	P_rad = 0.5*real(VandI(1)*conj(VandI(0)))/2+0.5*real(VandI(3)*conj(VandI(2)))/2+0.5*real(VandI(5)*conj(VandI(4)))/2;
	I.save("I_rv.txt",raw_ascii);
	V.save("V_r.txt",raw_ascii);*/
}

void rwg4_radiate_I()
{
	rwg3_main();
	cx_colvec Vv(EdgesTotal, 1), tempV, Ii;
	Vv.zeros();
	cx_mat Zz;
	rowvec fe = sort(FeedEdge, 1);
	if (FeedEdge.n_elem == 1) {
		Vv = -Z.col(FeedEdge(0));
		Vv.shed_row(FeedEdge(0));
		Zz = Z;
		Zz.shed_col(FeedEdge(0));
		Zz.shed_row(FeedEdge(0));
	}
	else if (FeedEdge.n_elem > 1) {
		for (int m = 0; m < FeedEdge.n_elem; m++)
			Vv = Vv - Z.col(FeedEdge(m));
		Zz = Z;
		tempV = V;
		for (int m = 0; m < fe.n_elem; m++) {
			Vv.shed_row(fe(m));
			Zz.shed_col(fe(m));
			Zz.shed_row(fe(m));
		}
	}
	cout << "\n" << "���������󷽳̣���ȴ�..." << endl;
	clock_t start = clock();
	Ii = solve(Zz, Vv);
	clock_t end = clock();
	cout << "\n" << "���󷽳����ʱ�䣺 " << (end - start) / 1000.0 << "s" << endl;
	for (int m = fe.n_elem - 1; m >= 0; m--) {
		cx_rowvec insert(1);
		insert(0) = complex<double>(1.0, 0.0);
		Ii.insert_rows(fe(m), insert);
	}
	Vv = tempV;

	V = Vv;
	I = Ii;
	cx_rowvec VandI(1, 2);// Vgap = 0, Igap = 0;
	VandI.zeros();

	for (int m = 0; m < FeedEdge.n_elem; m++) {
		VandI(0) = VandI(0) + as_scalar(Z.row(FeedEdge(m))*Ii);
		VandI(1) = VandI(1) + EdgeLength(FeedEdge(m));
	}
	P_rad = real(VandI(0)) / 2;
	Zin = VandI(0) / VandI(1) / VandI(1);

	cout << "\n" << "Zin: " << Zin << endl;
	//I.save("I_r.txt",raw_ascii);
	//V.save("V_r.txt",raw_ascii);
}

