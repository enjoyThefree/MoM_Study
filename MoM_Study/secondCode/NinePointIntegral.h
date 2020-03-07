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
extern vec EdgeLength;  //�����߳���
extern vec Vertex_pos;  // �������ε� �ǹ����ߵ� ���
extern vec Vertex_naga;

cube Center_;   //�ŵ���ֻ��ֺ��С����������

mat Pho_Pos;       //������element������ �ǹ����ߵ㵽���ĵ��������� 
mat Pho_Naga;      //���ڸ�element������ �ǹ����ߵ㵽���ĵ���������
cube Pho_Pos_;     //������element������ �ŵ㻮�ֺ� �ǹ����ߵ㵽 С���������ĵ���������
cube Pho_Naga_;    //���ڸ�element������ �ŵ㻮�ֺ� �ǹ����ߵ㵽 С���������ĵ���������

void NinepointIntegral() {

	InputData();
	cube center_(9, 3, elem_num);
	rowvec r1(1, 3), r2(1, 3), r3(1, 3);
	rowvec r12(1, 3), r23(1, 3), r13(1, 3);
	rowvec M(1, 3), c1(1, 3), c2(1, 3), c3(1, 3), c4(1, 3), c5(1, 3), c6(1, 3);

	int NoVertexP = 0; int NoVertexN = 0;   //�ǹ����ߵ� ���
	
	rowvec VertexP(1, 3);   // �������ε� �ǹ����ߵ� ������Ϣ
	rowvec VertexN(1, 3);   // �������ε� �ǹ����ߵ� ������Ϣ
	mat FreeVertex_Pmat(9, 3);    //�ŵ㻮�ֺ�,�������ε� ��չ�� �ǹ����ߵ� ������Ϣ
	mat FreeVertex_Nmat(9, 3);    //�ŵ㻮�ֺ�,�������ε� ��չ�� �ǹ����ߵ� ������Ϣ
	mat pho_pos(rwg_num, 3);     
	mat pho_nega(rwg_num, 3);    

	cube pho_pos_(9, 3, rwg_num);   
	cube pho_nega_(9, 3, rwg_num);

	for (int i = 0; i < elem_num; i++) {
		//���һ�� element ������������ꣻ
		r1 = element_list.slice(i).row(0);
		r2 = element_list.slice(i).row(1);
		r3 = element_list.slice(i).row(2);

		//�����������ߵ�����
		r12 = r2 - r1;
		r23 = r3 - r2;
		r13 = r3 - r1;

		// element ������
		M = Center.row(i);

		// element �����α��ϵ����ȷֵ�
		c1 = r1 + (1 / 3.0) * r12;
		c2 = r1 + (2 / 3.0) * r12;

		c3 = r2 + (1 / 3.0) * r23;
		c4 = r2 + (2 / 3.0) * r23;

		c5 = r1 + (1 / 3.0) * r13;
		c5 = r1 + (2 / 3.0) * r13;

		//����ÿ��С�����ε�����
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
		// ȡ�� ������Ϊ������rwg��element �ı��
		int NoPos = basis_list(i, 0);
		int NoNega = basis_list(i, 1);
		// ȡ�� �����������ε� �ǹ����ߵ� ���
		NoVertexP = Vertex_pos(i);
		NoVertexN = Vertex_naga(i);

		VertexP = element_list.slice(NoPos).row(NoVertexP);
		VertexN = element_list.slice(NoNega).row(NoVertexN);

		for (int i = 0; i < 9; i++) {
			FreeVertex_Pmat.row(i) = VertexP;
			FreeVertex_Nmat.row(i) = VertexN;
		}

		//�ǹ����ߵ㵽���ĵ��������� 
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