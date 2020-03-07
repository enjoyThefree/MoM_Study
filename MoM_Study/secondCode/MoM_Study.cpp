#include <iostream>
#include<armadillo>
#include"InputData.h"
#include"EFIE.h"
#include"NinePointIntegral.h"

using namespace arma;
using namespace std;
//using namespace mom;

mat basis_list;
cube element_list;
extern vec Area;     //单元面积
extern mat center;   //单元重心坐标
extern int elem_num;
extern int rwg_num;
//bool readData(mat& basis_list, mat& element_list);

const double pi = 3.1415926535898;
double mu = 4 * pi * 1e-7;                // 磁导率
double epsilon = 8.854187817 * 1e-12;     // 介电常数
double f = 0;                             // 定义入射波频率
double c = 1.0 / sqrt(mu * epsilon);     
double lamda = 0;                        // 波长
double omega = 0;                        // 角频率
double k = 0;                            // 波数                   
double N0 = 1.9e7;               //快速计算时格林函数采样点数
int FastGF = 0;                  //快速格林函数算法
double  N0 = 1.9e7;              //格林函数快速计算时采样点数
double  delta = 0;               //格林函数快速计算时采样间隔;
cx_rowvec Phase_base(N0 + 1);    //存放离散的相位

mat LoadMat(1, 3);
int main()
{
 
    //printf("Start readData\n");
    //readData(basis_list, element_list);
    ////InputData data_input;
    ////data_input.readData( basis_list, element_list);
    //basis_list.print("basis_list:");
    //element_list.print("element_list:");
    //printf("End readData\n");

    //printf("每个 element 面积\n");
    //area();
    //Area.print("Area:");
  
    //edgelength();
    //EdgeLength.print("edgeLength:");

    NinepointIntegral();
    Area.print("Area:");
    Center.print("Center:");
    Vertex_pos.print("Vertex_pos:");
    Pho_Pos.print("Pho_Pos:");
}

