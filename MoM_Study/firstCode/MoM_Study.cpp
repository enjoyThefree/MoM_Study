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

