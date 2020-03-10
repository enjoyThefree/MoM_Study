#pragma once
#include<iostream>
#include<armadillo>
#include<math.h>
#include<complex>
#include<time.h>
#include<Windows.h>
#include<cmath>
#include"NinePointIntegral.h"

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
extern vec EdgeLength;  //�����߳���
extern vec Vertex_pos;  // �������ε� �ǹ����ߵ� ���
extern vec Vertex_naga;

extern cube Center_;   //�ŵ���ֻ��ֺ��С����������

extern mat Pho_Pos;       //������element������ �ǹ����ߵ㵽���ĵ��������� 
extern mat Pho_Naga;      //���ڸ�element������ �ǹ����ߵ㵽���ĵ���������
extern cube Pho_Pos_;     //������element������ �ŵ㻮�ֺ� �ǹ����ߵ㵽 С���������ĵ���������
extern cube Pho_Naga_;    //���ڸ�element������ �ŵ㻮�ֺ� �ǹ����ߵ㵽 С���������ĵ���������

extern const double pi = 3.1415926535898;
extern double mu = 4 * pi * 1e-7;
extern double epsilon = 8.854187817 * 1e-12;
extern double f = 0;
extern double c = 1.0 / sqrt(mu * epsilon);
extern double lamda = 0;
extern double omega = 0;
extern double k = 0;
extern int FastGF ;    //���ٸ��ֺ����㷨
extern double delta;
extern cx_rowvec Phase_base;

extern mat LoadMat;
double constant1 = mu / (4 * pi);
double factor = 1 / 9.0;

cx_mat Z;

void FillmatZ()
{
	double FastPhi_index = 0;
	int Fastindex = 0;
	NinepointIntegral();

	complex<double> constant2(0, -1 / (4 * pi * omega * epsilon));
	cx_colvec Factor_A(rwg_num, 1), Factor_Fi(rwg_num, 1);

	cube Pho_P(9, 3, rwg_num), Pho_N(9, 3, rwg_num);        //   �ǹ����߶�Ӧ�㵽���� ��ʸ�� ��չ 9*3*rwg_num��ÿһ�㶼��ͬһ��rwg������һ��ԭʼ�����ε� pho��1*3 --> 9*3

	/*cx_mat z(rwg_num, rwg_num);*/
	rowvec Pos_index, Nega_index;   // ��RWG�� ����������element��� ����λ��
	cube RepmatCenter(9, 3, elem_num);   //element ������չ 9*3*elem_num
	cube RepmatCenter_(9, 3, elem_num);  //element�ŵ㻮�ֺ�����
	mat D(9, 3);   // elementԭ���������� ��ȥ �ŵ㻮�ֺ��С���������� ������
	mat R(elem_num, 9);     // R = |r - r'|  ; D�ı����� 
	cx_mat g(elem_num, 9);

	cx_mat gP(rwg_num, 9), gN(rwg_num, 9);  //���������ζ�Ӧ�ĸ��ֹ�ʽ

	cx_colvec Fi(rwg_num, 1), ZF(rwg_num, 1), Z1(rwg_num, 1);
	cx_colvec A(rwg_num, 1);
	cube R_PN(9, 3, rwg_num);     //������element������ �ŵ㻮�ֺ� �ǹ����ߵ㵽 С���������ĵ���������


	for (int i = 0; i < rwg_num; i++) {
		Factor_A(i) = complex<double>(0, factor * (omega * EdgeLength(i) / 4.0) * constant1);
		Factor_Fi(i) = complex<double>(factor * EdgeLength(i), 0) * constant2;

		mat temp_p(9, 3), temp_n(9, 3);  //�ݴ�ÿ��rwg��pho
		for (int j = 0; j < 9; j++) {
			temp_p.row(j) = Pho_Pos.row(i);
			temp_n.row(j) = Pho_Naga.row(i);
		}
		Pho_P.slice(i) = temp_p;     //��m��element��   �ǹ����߶�Ӧ�� �� ���� ��ʸ�� ��չ 9*3*rwg_num
		Pho_N.slice(i) = temp_n;
	}

	printf("\n ��ʼ��� Z����...\n");
	clock_t start = clock();

	cx_mat z(rwg_num, rwg_num);
	z = zeros<cx_mat>(rwg_num, rwg_num);

	//һ�� element������ Ҫ���������е� element������(�������Լ���������) ���и��ֺ�������
	for (int m = 0; m < elem_num; m++) {

		// ��RWG�� ����������element��� ����λ��
		int num_pos = 0;  	  // 
		int num_nega = 0;
		// �ж� ��(��m��)�������ǲ���rwg��������������
		for (int i = 0; i < rwg_num; i++) {
			if (basis_list(i, 0) == m) num_pos++;
			if (basis_list(i, 1) == m) num_nega++;
		}
		if (num_pos > 0) {
			int ip = 0;
			rowvec pos(num_pos);
			for (int i = 0; i < rwg_num; i++) {
				if (basis_list(i, 0) == m) {
					pos(ip) = i;
					ip++;
				}
			}
			Pos_index = pos;
		}
		if (num_nega > 0) {
			int in = 0;
			rowvec nega(num_nega);
			for (int i = 0; i < rwg_num; i++) {
				if (basis_list(i, 1) == m) {
					nega(in) = i;
					in++;
				}
			}
			Nega_index = nega;
		}

		/*������ֺ������� m ����� n ��element�ĸ��ֺ���*/
		//cube RepmatCenter(9, 3, rwg_num);   //element ������չ
		//cube RepmatCenter_(9, 3, rwg_num);  //element�ŵ㻮�ֺ�����
		for (int n = 0; n < elem_num; n++) {
			//element�ŵ㻮��������չ��һ��9��С������element���� --> 9*3*elem_num
			RepmatCenter_.slice(n) = Center_.slice(m);

			//���� n ��ԭʼelement��������չ�� 9*3*elem_num
			mat temp(9, 3);
			for (int i = 0; i < 9; i++) {
				temp.row(i) = Center.row(n);
			}
			RepmatCenter.slice(n) = temp;

			D = RepmatCenter.slice(n) - RepmatCenter_.slice(n);

			for (int i = 0; i < 9; i++) {
				R(n, i) = sqrt(as_scalar(D.row(i) * strans(D.row(i)))); //ÿһ��Ԫ�ض���һ��|r-r'|����
				if (FastGF == 0) {
					double kR = k * R(n, i);
					g(n, i) = exp(complex<double>(0, -kR)) / R(n, i);  //���ֺ���
				}
				else if (FastGF == 1) {
					double kR = k * R(n, i);
					FastPhi_index = (kR - 2 * pi * floor(kR / (2 * pi))) / delta;
					Fastindex = floor(FastPhi_index + 0.5);
					g(n, i) = Phase_base(Fastindex) / R(n, i);
				}

			}
		}

		for (int i = 0; i < rwg_num; i++) {
			gP.row(i) = g.row(basis_list(i, 0));
			gN.row(i) = g.row(basis_list(i, 1));
		}

		//sum( g(+) - g(-) )  
		Fi = sum(gP, 1) - sum(gN, 1);   //sum(A,dim),dim=0��ʾ����ӣ�dim=1��ʾ�����
		ZF = Factor_Fi % Fi; 

		int p_index = 0;
		if (num_pos > 0) {
			for (int i = 0; i < num_pos; i++) {
				p_index = Pos_index(i);
				for (int j = 0; j < rwg_num; j++) {
					R_PN.slice(j) = Pho_Pos_.slice(p_index);  
				}
				cube tempP_(9, 3, rwg_num), tempN_(9, 3, rwg_num);
				tempP_ = R_PN % Pho_P;    //
				tempN_ = R_PN % Pho_N;

				mat tempP(rwg_num, 9), tempN(rwg_num, 9);
				for (int j = 0; j < rwg_num; j++) {
					tempP.row(j) = strans(sum(tempP_.slice(j), 1));   // ��һ��element�� 9��С�����ε� ����������ӣ�
					tempN.row(j) = strans(sum(tempN_.slice(j), 1));

				}
				A = sum(gP % tempP, 1) + sum(gN % tempN, 1);
				Z1 = Factor_A % A;
				z.col(p_index) = z.col(p_index) + EdgeLength(p_index) * (Z1 + ZF);
			}
		}

		int n_index = 0;
		if (num_nega > 0) {
			for (int i = 0; i < num_nega; i++) {
				n_index = Nega_index(i);
				for (int j = 0; j < rwg_num; j++) {
					R_PN.slice(j) = Pho_Pos_.slice(n_index);
				}
				cube tempP_(9, 3, rwg_num), tempN_(9, 3, rwg_num);
				tempP_ = R_PN % Pho_P;    //
				tempN_ = R_PN % Pho_N;

				mat tempP(rwg_num, 9), tempN(rwg_num, 9);
				for (int j = 0; j < rwg_num; j++) {
					tempP.row(j) = strans(sum(tempP_.slice(j), 1));   // ��һ��element�� 9��С�����ε� ����������ӣ�
					tempN.row(j) = strans(sum(tempN_.slice(j), 1));

				}
				A = sum(gP % tempP, 1) + sum(gN % tempN, 1);
				Z1 = Factor_A % A;
				z.col(n_index) = z.col(n_index) + EdgeLength(n_index) * (Z1 + ZF);
			}
		}

	}

	for (int i = 0; i < rwg_num; i++) {
		int element_index = basis_list(i, 0);  //ȡ�� element���
		rowvec temp_node = CommonEdgeP.row(i);  //rwg�����ߵ�������
		// �����ߵ��е�����
		rowvec temp_center = (element_list.slice(element_index).row(temp_node(0)) + element_list.slice(element_index).row(temp_node(1))) / 2.0;

		if (norm(temp_center - LoadMat) < 1e-8) {
			z(i, i) = z(i, i) + pow(EdgeLength(i), 2) * 50;  //50: complex<double>(0,omega*50.0)
		}
	}

	clock_t end = clock();
	printf("\n Z�������ʱ��Ϊ \t");
	cout << (end - start) << "s" << endl;

	Z = z;
} 