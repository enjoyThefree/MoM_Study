#include <iostream>
#include<armadillo>
#include"InputData.h"
#include"EFIE.h"
#include"NinePointIntegral.h"
#include"FillZmat.h"
#include"SolveMat.h"

using namespace arma;
using namespace std;
//using namespace mom;
void storeMoRCS();
cx_rowvec CrossComplex(cx_rowvec a, cx_rowvec b);
void DipoleSFarField(rowvec Point);
float MonoRCS(rowvec PRE_ObservationPoint, rowvec pol, double FarR);
void DipoleModeling();

mat basis_list;
cube element_list;
extern vec Area;     //单元面积
extern mat Center;   //单元重心坐标
extern int elem_num;
extern int rwg_num;
extern cube Center_;   //九点积分划分后的小三角形重心

extern mat Pho_Pos;       //对于正element三角形 非公共边点到重心的向量坐标 
extern mat Pho_Naga;      //对于负element三角形 非公共边点到重心的向量坐标
extern cube Pho_Pos_;     //对于正element三角形 九点划分后 非公共边点到 小三角形重心的向量坐标
extern cube Pho_Naga_;    //对于负element三角形 九点划分后 非公共边点到 小三角形重心的向量坐标

extern cx_mat Z;
extern cx_colvec I, V;

mat DipoleCenter;
cx_mat DipoleMoment;
cx_rowvec HField(1, 3), EField(1, 3);
rowvec BIRCS, MoRCS, RadiRCS;
//bool readData(mat& basis_list, mat& element_list);

/* ****************  计算参数设置  ******************************/
const double pi = 3.1415926535898;
double mu = 4 * pi * 1e-7;                // 磁导率
double epsilon = 8.854187817 * 1e-12;     // 介电常数
double f = 0;                             // 定义入射波频率
double c = 1.0 / sqrt(mu * epsilon);     
double eta = sqrt(mu / epsilon);         // 波阻抗
double lamda = 0;                        // 波长
double omega = 0;                        // 角频率
double k = 0;                            // 波数                   

int FastGF = 0;                  //快速格林函数算法 是否采用格林函数快速计算方法：1-是  0-否
double  N0 = 1.9e7;              //格林函数快速计算时采样点数
double  delta = 0;               //格林函数快速计算时采样间隔;
cx_rowvec Phase_base(N0 + 1);    //存放离散的相位

rowvec Pol(1, 3), dz(1, 3);      //定义入射波的极化方向和入射方向(only for r_s = 1)
double  Eamp = 0;                //定义入射波电场幅度
rowvec  FeedPoint(1, 3);         //定义电压源的馈电点(only for r_s = 2)
rowvec  FeedPlus(1, 3);          //定义馈电点的正极(only for r_s = 2)
rowvec FeedEdge;                 //定义馈电边的EdgeElement编号
double FeedVoltage = 0.0;        //定义馈电电压

int     r_s = 0;                 //定义求解问题的类型：1-scattering  2-radiating
int     RCS_type = 1;            //RCS类型：1-单站  2-双站
int     RCS_plane = 0;            //待计算的RCS截面：1-xoz  2-yoz  3-xoy

float Emag = 0;                   //定义入射波电场幅度
double Distance = 0;              //计算距离

mat LoadMat(1, 3);

int main()
{
    delta = 2 * pi / N0;
    int phase_count = 0;
    for (double n = 0; n < N0; n++) {
        double phi = delta * n;
        Phase_base(phase_count) = exp(complex<double>(0, -phi));
        phase_count++;
    }
    Phase_base(N0) = Phase_base(0);

    /*------------------------------------*/
    f = 1.0e9;   //入射波频率
    r_s = 1;     
    lamda = c / f;
    omega = 2.0 * pi * f;
    k = omega / c;

    if (r_s == 1) {   // 散射问题
        Emag = 1;
       // RCS_type = 1;
        if (RCS_type == 1) {
            Distance = 1000;
            float cal_angle = 0;
            rowvec morcs(1, 361);
            FillmatZ();
            /* 求阻抗矩阵的逆 */
            cx_mat Z_inv = inv(Z);

            float theta = 90 * pi / 180;
            int count = 0;

            clock_t start = clock();
            for (float phi = -180; phi <= 180; phi++) {
                cal_angle = phi * pi / 180;   //将角度用 pi 表示
                dz << -sin(theta) * cos(cal_angle) << -sin(theta) * sin(cal_angle) << -cos(theta) << endr;
                Pol << 0.0 << 0.0 << 1.0 << endr;
                Pol = Pol * Emag;
                SolveMat_V();
                I = Z_inv * V;
                DipoleModeling();
                morcs(count) = MonoRCS(dz, Pol, Distance);
                count++;
            }
            MoRCS = morcs;
            clock_t end = clock();

            cout << endl << "计算时间为：" << (end - start) / 1000 << "s" << endl;
            cout << endl << "正在保存RCS数据 。。。" << endl;
            void storeMoRCS();
        }
    }


   // FillmatZ();
   //// NinepointIntegral();
   // Area.print("Area:");
   // Center.print("Center:");
   // Vertex_pos.print("Vertex_pos:");
   // Pho_Pos.print("Pho_Pos:");

   // Z.print("Z:");
}

/* 将基函数仿真成 偶极子模型 */
void DipoleModeling() {
    rowvec Point1, Point2;
    Point1 = zeros<rowvec>(1, 3);
    Point2 = zeros<rowvec>(1, 3);
    mat dipoleCenter(rwg_num, 3);
    cx_mat dipoleMoment(rwg_num, 3);

    for (int i = 0; i < rwg_num; i++) {
        Point1 = Center.row(basis_list(i, 0));
        Point2 = Center.row(basis_list(i, 1));

        dipoleCenter.row(i) = 0.5 * (Point1 + Point2);
        dipoleMoment.row(i) = EdgeLength(i) * I(i) * (-1.0 * Point1 + Point2);

    }
    DipoleCenter = dipoleCenter;
    DipoleMoment = dipoleMoment;
}

/* 计算单站RCS */
float MonoRCS(rowvec PRE_ObservationPoint, rowvec pol, double FarR) {
    double EIabs = 0, ConstantRCS = 0;
    float monorcs = 0; 
    rowvec ObservationPoint = zeros<rowvec>(1.3);
    ObservationPoint = -1.0 * FarR * PRE_ObservationPoint;

    EIabs = sqrt(as_scalar(pol * strans(pol)));
    ConstantRCS = 4 * pi * FarR * FarR / (EIabs * EIabs);
    DipoleSFarField(ObservationPoint);
    double ESabsPower2 = 0.0;
    for (int i = 0; i < 3; i++)
        ESabsPower2 = ESabsPower2 + abs(EField(i)) * abs(EField(i));
    monorcs = 10 * log10(ConstantRCS * ESabsPower2);
    return monorcs;
}

/* 计算所有电偶极子在空间某一点叠加的远场 */
void DipoleSFarField(rowvec Point) {
    double ConstantE = eta / (4.0 * pi);
    complex<double>ConstantH(0, k / (4.0 * pi));
    complex<double>K2(0, k);

    HField = zeros<cx_rowvec>(1, 3);
    EField = zeros<cx_rowvec>(1, 3);
    
    rowvec Center1(1, 3), R1(1, 3);
    cx_rowvec Moment(1, 3), M(1, 3), MomentCrossR(1, 3), R1Complex(1, 3);
    double PointRM, PointRM2;
    complex<double> EXP1(0, 0);
    complex<double> C2(0, 0);

    for (int i = 0; i < rwg_num; i++) {
        Moment = DipoleMoment.row(i);
        Center1 = DipoleCenter.row(i);
        R1 = Point - Center1;
        R1Complex = R1 * (complex<double>(1.0, 0.0));
        PointRM2 = as_scalar(R1 * strans(R1));
        PointRM = sqrt(PointRM2);
        EXP1 = exp(-K2 * PointRM);
        C2 = 1.0 / PointRM2 * (1.0 + 1.0 / (K2 * PointRM));
        M = as_scalar(R1 * strans(Moment)) * R1 / PointRM2;
        MomentCrossR = CrossComplex(Moment, R1Complex);
        HField = HField + ConstantH * MomentCrossR * C2 * EXP1;
        EField = EField + ConstantE * ((M - Moment) * (K2 / PointRM + C2) + 2 * M * C2) * EXP1;
    }
}

/* 计算复向量的叉乘 */
cx_rowvec CrossComplex(cx_rowvec a, cx_rowvec b) {
    cx_rowvec c(1, 3);
    c(0) = a(1) * b(2) - a(2) * b(1);
    c(1) = a(2) * b(0) - a(0) * b(2);
    c(2) = a(0) * b(1) - a(1) * b(0);
    return c;
}

/* 保存数据到文件 */
void storeMoRCS()
{
    string filename = "MoRCS.txt";
    string line;
    int i;

    ofstream fout(filename.c_str());

    if (!fout.is_open()) {
        cout << "打开文件失败" << endl;
    }
    
    for (i = 0; i < MoRCS.size(); i++) {
        fout << MoRCS[i] << " ";
    }
    fout.close();
}