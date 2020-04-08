#pragma once
#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<armadillo>
#include <sstream>
#include<math.h>

using namespace std;
using namespace arma;

extern mat basis_list;
extern cube element_list;
int elem_num;
int rwg_num;
mat Center;
vec Area;

mat CommonEdgeP;  //������RWG ��ɹ����ߵ���������Ϣ
mat CommonEdgeN;  //������RWG ��ɹ����ߵ���������Ϣ
vec EdgeLength;  //�����߳���

vec Vertex_pos;  // �������ε� �ǹ����ߵ� ���
vec Vertex_naga; // �������ε� �ǹ����ߵ� ���

//mat basis_list;
//mat element_list;

namespace mom 
{
	class InputData {
	public:
		InputData(){}
		//bool readData(string &fileName,mat &basis_list,mat &element_list);
		bool readData(mat &basis_list, cube &element_list);
	};

}

bool InputStream(ifstream& in, stringstream& ss) {
	string str;
	ss.clear();
	getline(in, str);
	ss << str.substr(0, str.find("#"));
	return true;
}

bool readData(mat& basis_list, cube& element_list) {
	string ignore;
	stringstream ss;
	/*��ȡrwg:�����Լ����ڵ�����element*/
	int num = 0, basis_p = 0, basis_n = 0;
	ifstream basis_in("rwg.txt", ios::in);
	if (!basis_in.is_open())
		cout << "open rwg.txt file failure";

	//basis_in >> num;
	InputStream(basis_in, ss);
	ss >> num;
	rwg_num = num;
	basis_list.set_size(num, 2);
	for (int k = 0; k < num; k++) {
		InputStream(basis_in, ss);
		ss >> basis_p >> basis_n >> ignore >> ignore;
		basis_list(k, 0) = basis_p;
		basis_list(k, 1) = basis_n;
	}
	//int k;
	//while (num--) {
	//	//basis_in >> basis_p >> basis_n >> ignore >> ignore;
	//	InputStream(basis_in, ss);
	//	ss >> basis_p >> basis_n >> ignore >> ignore;
	//	basis_list(k, 0) = basis_p;
	//	basis_list(k, 1) = basis_n;
	//	k++;
	//}
	basis_in.close();

	/*��ȡ����,element���*/
	ifstream element_in("tri.txt", ios::in);
	double x, y, z;
	//mat vertex_list(1, 3);
	if (!element_in.is_open())
		cout << "open tri.txt file failure";
	InputStream(element_in, ss);
	ss >> num;
	elem_num = num;
	element_list.set_size(3, 3, num);  //��һ��Ϊ��һ��element
	for (int i = 0; i < num; i++) {
		mat tem(3, 3);  //��һ�У�p1�� ...  ����
		for (int j = 0; j < 3; j++) {
			InputStream(element_in, ss);
			ss >> x >> y >> z;
			tem(j, 0) = x;
			tem(j, 1) = y;
			tem(j, 2) = z;
		}
		element_list.slice(i) = tem;
	}
	element_in.close();

	return true;
}

/*�����������*/
double* cross(double x1[3], double x2[3])
{
	double* x3 = new double[3];
	x3[0] = x1[1] * x2[2] - x1[2] * x2[1];
	x3[1] = x1[2] * x2[0] - x1[0] * x2[2];
	x3[2] = x1[0] * x2[1] - x1[1] * x2[0];
	return x3;
}
/*����ÿ������Ƭ���������ص�*/
void area()
{
	vec area(elem_num);
	mat center(elem_num, 3);
	double vector1[3], vector2[3];   //������
	double * vector3 = new double[3]; //��˽��
	for (int i = 0; i < elem_num; i++) {
		for (int j = 0; j < 3; j++) {
			vector1[j] = element_list.slice(i)(0, j) - element_list.slice(i)(1, j);
			vector2[j] = element_list.slice(i)(2, j) - element_list.slice(i)(1, j);
		}
		vector3 = cross(vector1, vector2); 
		area(i) = 0.5 * sqrt(pow(vector3[0], 2) + pow(vector3[1], 2) + pow(vector3[2], 2));
		/*center(i, 0) = (1.0 / 3.0) * (element_list.slice(i)(0, 0) + element_list.slice(i)(1, 0) + element_list.slice(i)(2, 0));
		center(i, 1) = (1.0 / 3.0) * (element_list.slice(i)(0, 1) + element_list.slice(i)(1, 1) + element_list.slice(i)(2, 1));
		center(i, 2) = (1.0 / 3.0) * (element_list.slice(i)(0, 2) + element_list.slice(i)(1, 2) + element_list.slice(i)(2, 2));*/
		for (size_t k = 0; k < 3; k++) {
			center(i, k) = (1.0 / 3.0) * (element_list.slice(i)(0, k) + element_list.slice(i)(1, k) + element_list.slice(i)(2, k));
		}
	}
	delete[] vector3;
	Area = area;
	Center = center;
}

/*Ѱ����ɹ����ߵ������� �Լ� �ǹ����ߵĵ�*/
void commonEdge() {
	mat commonEdge_p(rwg_num, 2);
	mat commonEdge_n(rwg_num, 2);
	vec vertex_pos(rwg_num);
	vec vertex_nega(rwg_num);
	for (int i = 0; i < rwg_num; i++) {
		int t = 0;
		for (int m = 0; m < 3; m++) {
			for (int n = 0; n < 3; n++) {
				//slice ѡ�㣬��ѡ element    (basis_list ����������˵����� element)
			// (m,0)  m:���element�ĵ�m����   0:x����   1��y����  2��z����
				if (element_list.slice(basis_list(i, 0))(m, 0) == element_list.slice(basis_list(i, 1))(n, 0) &&
					element_list.slice(basis_list(i, 0))(m, 1) == element_list.slice(basis_list(i, 1))(n, 1) &&
					element_list.slice(basis_list(i, 0))(m, 2) == element_list.slice(basis_list(i, 1))(n, 2))
				{
					commonEdge_p(i, t) = m;
					commonEdge_n(i, t) = n;
					t++;
				}
				//if (t > 2) break;
			}
			if (t > 2) {
				printf("��ͬ��������2���˳�");
				break;
			}
		}
		for (int k = 0; k < 3; k++) {
			if (k != commonEdge_p(i,0) && k != commonEdge_p(i,1)) {
				vertex_pos(i) = k;
			}
			if (k != commonEdge_n(i, 0) && k != commonEdge_n(i, 1)) {
				vertex_nega(i) = k;
			}
		}

		//if (t > 2) break;
	}
	CommonEdgeP = commonEdge_p;
	CommonEdgeN = commonEdge_n;
	Vertex_pos = vertex_pos;
	Vertex_naga = vertex_nega;
}
void edgelength()
{
	commonEdge();
	vec edgeLength(rwg_num);
	rowvec vector1(1, 3), vector2(1, 3), vector3(1, 3);
	for (int i = 0; i < rwg_num; i++) {
		// ȡ�� ������Ϊ����element �ı��
		int NoPos = basis_list(i, 0);
		vector1 = element_list.slice(NoPos).row(CommonEdgeP(i, 0)) - element_list.slice(NoPos).row(CommonEdgeP(i, 1));
		edgeLength(i) = sqrt(as_scalar(vector1 * strans(vector1)));
	}
	EdgeLength = edgeLength;
}

void InputData() {
	readData(basis_list, element_list);
	area();
	edgelength();
}