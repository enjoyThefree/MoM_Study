/*����õ�Ŀ��������ϵ��֮�󣬾Ϳ��Ը������ɿռ����Դ��ɢ�䳡��ʽ
����Ŀ�������䲨�����µĿռ�ɢ�䳡�������������*/

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <complex>
#include <math.h>
#include <time.h>
#include <windows.h>
#include <cmath>
//#include "RCS.h"
#include "rwg1.h"
#include "rwg2.h"
#include "rwg3.h"
#include "rwg4.h"
#include "armadillo"

using namespace arma;
using namespace std;
const double pi = 3.14151954;

extern float  P_rad;
extern int    TrianglesTotal, EdgesTotal;             //�ֱ��ʾ������Ŀ���������Ƭ���Լ�RWG��������
extern mat    t_data;                                 //�þ���洢Triangle�ڵ�����Ϣ
extern mat    p_data;                                 //�þ���洢Triangle�ڵ�������Ϣ
extern mat    EdgeElement;                            //�þ���洢EdgeElement�ڵ�����Ϣ
extern vec    TrianglePlus;                           //�þ���洢EdgeElement��Ӧ���������α����Ϣ
extern vec    TriangleMinus;                          //�þ���洢EdgeElement��Ӧ�ĸ������α����Ϣ                    
extern vec    EdgeLength;                             //�þ���洢EdgeElement�ĳ���
extern mat    Center;                                 //�þ���洢Triangle��Ԫ����������
extern vec    Area;                                   //������洢Triangle��Ԫ�����
extern mat    Pho_Plus;                               //�þ���洢EdgeElement��Ӧ����������RWG
extern mat    Pho_Minus;                              //�þ���洢EdgeElement��Ӧ�ĸ�������RWG
extern cube   Center_;                                //�þ���洢Triangle��Ԫ�����ŵ㻮�ֺ���������ε���������
extern cube   Pho_Plus_;                              //�þ���洢EdgeElement��Ӧ���������ξ����ŵ㻮�ֺ��RWG
extern cube   Pho_Minus_;                             //�þ���洢EdgeElement��Ӧ�ĸ������ξ����ŵ㻮�ֺ��RWG
extern float  Z_fill_time;

extern cx_mat    Z;
extern cx_colvec I, V;
extern complex<double> Zin;

mat        DipoleCenter;
cx_mat     DipoleMoment;
cx_rowvec  HField(1, 3), EField(1, 3);
rowvec     BiRCS, MoRCS, RadiRCS;

/*----------------------�����������----------------------------*/
/*--------------------------------------------------------------*/
string  t_filename = "microstrip_t.txt";
string  p_filename = "microstrip_p.txt";

double mu = 4 * pi * 1e-7;                // �ŵ���
double epsilon = 8.854187817 * 1e-12;     // ��糣��

int     r_s = 0;                                      //���������������ͣ�1-scattering  2-radiating
int     RCS_type = 0;                                 //RCS���ͣ�1-��վ  2-˫վ
int     RCS_plane = 0;                                //�������RCS���棺1-xoz  2-yoz  3-xoy
int     SimpleGF = 0;                                 //�Ƿ���ø��ֺ������ټ��㷽����1-��  0-��
double  N0 = 1.9e7;                                   //���ֺ������ټ���ʱ��������
double  delta = 0;                                    //���ֺ������ټ���ʱ�������;
double  Distance = 0;                                 //�������
float   f = 0;                                        //�������䲨Ƶ��
rowvec  Pol(1, 3), dz(1, 3);                             //�������䲨�ļ�����������䷽��(only for r_s = 1)
float   Emag = 0;                                     //�������䲨�糡����
rowvec  FeedPoint(1, 3);                               //�����ѹԴ�������(only for r_s = 2)
rowvec  FeedPlus(1, 3);                                //��������������(only for r_s = 2)
rowvec  FeedEdge;                                     //��������ߵ�EdgeElement���
float   FeedVoltage = 0.0;                            //���������ѹ
float   c = 1.0 / sqrt(epsilon*mu);
float   eta = sqrt(mu / epsilon);
float   lamda = 0;
float   omega = 0;
double  k = 0;
cx_rowvec g_base(N0 + 1);
/*--------------------------------------------------------------*/
mat FeedMat(1, 3), FeedPlusMat(1, 3), LoadMat(1, 3);
void main()
{
	delta = 2 * pi / N0;
	//---------build the value list for phase item-----------//
	int g_count = 0;
	for (double m = 0; m < N0; m++)
	{
		double phi = delta * m;
		g_base(g_count) = exp(complex<double>(0, -phi));
		g_count++;
	}
	g_base(N0) = g_base(0);
	//-------------------------------------------------------//
	double ElapsedTime = 0;

	f = 1.0e9;
	r_s = 1;
	lamda = c / f;
	omega = 2.0*pi*f;
	k = omega / c;

	if (r_s == 1) {
		Emag = 1;
		//�������ɢ������,���ô������RCS����
		RCS_type = 2;
		if (RCS_type == 1) {
			Distance = 1000;
			float cal_angle = 0;
			rowvec morcs(1, 361);
			rwg3_main();
			cout << "\n" << "�������迹������棬��ȴ�..." << endl;
			clock_t start1 = clock();
			cx_mat Z_inv = inv(Z);
			clock_t end1 = clock();
			cout << "\n" << "�迹��������ʱ�䣺 " << (end1 - start1) / 1000.0 << "s" << endl;

			float theta = (90)*pi / 180;
			int count = 0;
			clock_t start = clock();
			for (float phi = -180; phi <= 180; phi++)
			{
				cal_angle = phi * pi / 180;
				dz << -sin(theta)*cos(cal_angle) << -sin(theta)*sin(cal_angle) << -cos(theta) << endr;
				Pol << 0.0 << 0.0 << 1.0 << endr;
				Pol = Pol * Emag;
				rwg4_monoV();
				I = Z_inv * V;
				DipoleModeling();
				morcs(count) = MonoRCS(dz, Pol, Distance);
				count++;
			}
			MoRCS = morcs;
			clock_t end = clock();
			ElapsedTime = end - start;
			cout << "\n" << "�ܼ���ʱ�䣺" << ElapsedTime / 1000.0 << "s" << endl;
			cout << "\n" << "�������RCS���ݣ���ȴ�..." << "\n" << endl;
			MoRCS.save("MoRCS.txt", raw_ascii);
		}
		else if (RCS_type == 2) {
			//���ù۲��漰�������
			RCS_plane = 1;
			Distance = 1000;
			//�������䲨����ͼ���
			dz << 0.0 << 0.0 << -1.0 << endr;
			Pol << 1.0 << 0.0 << 0.0 << endr;
			Pol = Pol * Emag;
			/*---------------����˫վRCS�Ĺ۲���------------------------*/
			int AngleNum = 361;
			rowvec AngleValue(1, AngleNum);
			rowvec Theta, Phi;
			for (int ia = 0; ia < AngleNum; ia++)
				AngleValue(ia) = (ia - 180)*pi / 180.0;
			if (RCS_plane == 1) {
				rowvec theta(AngleNum), phi(1);
				phi(0) = 0.0*pi / 180.0;
				theta = AngleValue;
				Theta = theta;
				Phi = phi;
			}
			else if (RCS_plane == 2) {
				rowvec theta(AngleNum), phi(1);
				phi(0) = 90.0*pi / 180.0;
				theta = AngleValue;
				Theta = theta;
				Phi = phi;
			}
			else if (RCS_plane == 3) {
				rowvec theta(1), phi(AngleNum);
				theta(0) = 90.0*pi / 180.0;
				phi = AngleValue;
				Theta = theta;
				Phi = phi;
			}
			/*--------------------------------------------------------------*/
			clock_t start = clock();
			rwg4_main();
			DipoleModeling();
			BioRCS(AngleNum, Theta, Phi, Distance, RCS_plane);
			clock_t end = clock();
			ElapsedTime = end - start;
			cout << "\n" << "�ܼ���ʱ�䣺" << ElapsedTime / 1000.0 << "s" << endl;
			cout << "\n" << "�������RCS���ݣ���ȴ�..." << "\n" << endl;
			BiRCS.save("BiRCS.txt", raw_ascii);
			I.save("I.txt", raw_ascii);
		}
	}
	else if (r_s == 2) {
		FeedMat << 0.0 << 0.0 << 0.0 << endr;
		FeedPlusMat << 0.0 << 0.02 << 0.25 << endr;
		LoadMat << 0.0 << 0.0 << 0.5 << endr;
		FeedVoltage = 1.0;
		//���ù۲��漰�������
		RCS_plane = 1;
		Distance = 1000;
		int AngleNum = 361;
		rowvec AngleValue(1, AngleNum);
		rowvec Theta, Phi;
		for (int ia = 0; ia < AngleNum; ia++)
			AngleValue(ia) = (ia - 180)*pi / 180.0;
		if (RCS_plane == 1) {
			rowvec theta(AngleNum), phi(1);
			phi(0) = 0.0*pi / 180.0;
			theta = AngleValue;
			Theta = theta;
			Phi = phi;
		}
		else if (RCS_plane == 2) {
			rowvec theta(AngleNum), phi(1);
			phi(0) = 90.0*pi / 180.0;
			theta = AngleValue;
			Theta = theta;
			Phi = phi;
		}
		else if (RCS_plane == 3) {
			rowvec theta(1), phi(AngleNum);
			theta(0) = 90.0*pi / 180.0;
			phi = AngleValue;
			Theta = theta;
			Phi = phi;
		}
		clock_t start = clock();
		rwg4_radiate_V();
		DipoleModeling();
		RadiateGain(AngleNum, Theta, Phi, Distance, RCS_plane);
		clock_t end = clock();
		ElapsedTime = end - start;
		cout << "\n" << "�ܼ���ʱ�䣺" << ElapsedTime / 1000.0 << "s" << endl;
		//cout << "\n" << "�������RCS���ݣ���ȴ�..." << "\n" << endl;
		RadiRCS.save("RadiRCS.txt", raw_ascii);
	}
}
/***********************************************************/
/*              Function��DipoleModeling()                 */
/***********************************************************/
/* 1.������������ɵ�ż����ģ��                            */
/***********************************************************/
void DipoleModeling()
{
	rowvec Point1, Point2;
	Point1 = zeros<rowvec>(1, 3);
	Point2 = zeros<rowvec>(1, 3);
	mat dipoleCenter(EdgesTotal, 3);
	cx_mat dipoleMoment(EdgesTotal, 3);

	for (int ia = 0; ia < EdgesTotal; ia++)
	{
		Point1 = Center.row(TrianglePlus(ia));
		Point2 = Center.row(TriangleMinus(ia));
		dipoleCenter.row(ia) = 0.5*(Point1 + Point2);
		dipoleMoment.row(ia) = EdgeLength(ia)*I(ia)*((-1.0)*Point1 + Point2);
	}
	DipoleCenter = dipoleCenter;
	DipoleMoment = dipoleMoment;
}
/**************************************************************************/
/*       Function��cx_rowvec CrossComplex(cx_rowvec a, cx_rowvec b)       */
/**************************************************************************/
/* 1.���㸴�����Ĳ��                                                     */
/**************************************************************************/
cx_rowvec CrossComplex(cx_rowvec a, cx_rowvec b)
{
	cx_rowvec c(1, 3);
	c(0) = a(1)*b(2) - a(2)*b(1);
	c(1) = a(2)*b(0) - a(0)*b(2);
	c(2) = a(0)*b(1) - a(1)*b(0);
	return c;
}
/**************************************************************************/
/*               Function��void DipoleSFarField(rowvec Point)             */
/**************************************************************************/
/* 1.�������е�ż�����ڿռ�ĳһ��������ɵ�Զ��                           */
/**************************************************************************/
void DipoleSFarField(rowvec Point)
{
	double ConstantE = eta / (4.0*pi);
	complex<double> ConstantH(0, k / (4 * pi));
	complex<double> K2(0, k);

	HField = zeros<cx_rowvec>(1, 3);
	EField = zeros<cx_rowvec>(1, 3);

	rowvec Center1(1, 3), R1(1, 3);
	cx_rowvec Moment(1, 3), M(1, 3), MomentCrossR(1, 3), R1Complex(1, 3);
	double PointRM, PointRM2;
	complex<double> EXP1(0, 0);
	complex<double> C2(0, 0);

	for (int ia = 0; ia < EdgesTotal; ia++)
	{
		Moment = DipoleMoment.row(ia);
		Center1 = DipoleCenter.row(ia);
		R1 = Point - Center1;
		R1Complex = R1 * (complex<double>(1.0, 0.0));
		PointRM2 = as_scalar(R1*strans(R1));
		PointRM = sqrt(PointRM2);
		EXP1 = exp(-K2 * PointRM);
		C2 = 1.0 / PointRM2 * (1.0 + 1.0 / (K2*PointRM));
		M = as_scalar(R1*strans(Moment))*R1 / PointRM2;
		MomentCrossR = CrossComplex(Moment, R1Complex);
		HField = HField + ConstantH * MomentCrossR*C2*EXP1;
		EField = EField + ConstantE * ((M - Moment)*(K2 / PointRM + C2) + 2 * M*C2)*EXP1;
	}
}
/**************************************************************************/
/*   void BioRCS(int angleNum,rowvec angleValue,double FarR,int VV_or_HH) */
/**************************************************************************/
/* 1.����˫վRCS                                                          */
/**************************************************************************/
void BioRCS(int angleNum, rowvec theta, rowvec phi, double FarR, int RcsType)
{
	double EIabs = 0, ConstantRCS = 0;
	double xory = 0, z0 = 0;
	rowvec biRCS(angleNum);
	rowvec ObservationPoint = zeros<rowvec>(1, 3);

	EIabs = sqrt(as_scalar(Pol*strans(Pol)));
	//ConstantRCS = 4*pi*FarR*FarR/(EIabs*EIabs*lamda*lamda);   //���ּ��㷽ʽΪ�Բ�����һ�����RCS
	ConstantRCS = 4 * pi*FarR*FarR / (EIabs*EIabs);
	for (int ia = 0; ia < theta.n_elem; ia++)
	{
		xory = FarR * sin(theta(ia));
		z0 = FarR * cos(theta(ia));
		ObservationPoint(2) = z0;
		if (RcsType == 1) //xoz
		{
			ObservationPoint(0) = xory;
			ObservationPoint(1) = 0.0;
			DipoleSFarField(ObservationPoint);

			double ESabsPower2 = 0.0;
			for (int i = 0; i < 3; i++)
			{
				ESabsPower2 = ESabsPower2 + abs(EField(i))*abs(EField(i));
			}
			biRCS(ia) = 10 * log10(ConstantRCS*ESabsPower2);
		}
		else if (RcsType == 2) //yoz
		{
			ObservationPoint(0) = 0.0;
			ObservationPoint(1) = xory;
			DipoleSFarField(ObservationPoint);

			double ESabsPower2 = 0.0;
			for (int i = 0; i < 3; i++)
			{
				ESabsPower2 = ESabsPower2 + abs(EField(i))*abs(EField(i));
			}
			biRCS(ia) = 10 * log10(ConstantRCS*ESabsPower2);
		}
		else if (RcsType == 3) //xoz 
		{
			for (int ik = 0; ik < phi.n_elem; ik++)
			{
				ObservationPoint(0) = xory * cos(phi(ik));
				ObservationPoint(1) = xory * sin(phi(ik));
				DipoleSFarField(ObservationPoint);
				double ESabsPower2 = 0.0;
				for (int i = 0; i < 3; i++)
				{
					ESabsPower2 = ESabsPower2 + abs(EField(i))*abs(EField(i));
				}
				biRCS(ik) = 10 * log10(ConstantRCS*ESabsPower2);
			}

		}
	}
	BiRCS = biRCS;
}
/**************************************************************************/
/*    float MonoRCS(rowvec PRE_ObservationPoint,rowvec pol,double FarR)   */
/**************************************************************************/
/* 1.���㵥վRCS                                                          */
/**************************************************************************/
float MonoRCS(rowvec PRE_ObservationPoint, rowvec pol, double FarR)
{
	double EIabs = 0, ConstantRCS = 0;
	float monorcs = 0;
	rowvec ObservationPoint = zeros<rowvec>(1, 3);
	ObservationPoint = -1.0*FarR*PRE_ObservationPoint;

	EIabs = sqrt(as_scalar(pol*strans(pol)));
	//ConstantRCS = 4*pi*FarR*FarR/(EIabs*EIabs*lamda*lamda);   //���ּ��㷽ʽΪ�Բ�����һ�����RCS
	ConstantRCS = 4 * pi*FarR*FarR / (EIabs*EIabs);
	DipoleSFarField(ObservationPoint);
	double ESabsPower2 = 0.0;
	for (int i = 0; i < 3; i++)
		ESabsPower2 = ESabsPower2 + abs(EField(i))*abs(EField(i));
	monorcs = 10 * log10(ConstantRCS*ESabsPower2);
	return monorcs;
}
/*****************************************************************************/
/* void RadiateGain(int angleNum,rowvec angleValue,double FarR,int VV_or_HH) */
/*****************************************************************************/
/* 1.������䳡RCS                                                           */
/*****************************************************************************/
void RadiateGain(int angleNum, rowvec theta, rowvec phi, double FarR, int RcsType)
{
	double EIabs = 0, ConstantRCS = 0;
	double xory = 0, z0 = 0;
	rowvec radiRCS(angleNum);
	rowvec ObservationPoint = zeros<rowvec>(1, 3);

	for (int ia = 0; ia < theta.n_elem; ia++)
	{
		xory = FarR * sin(theta(ia));
		z0 = FarR * cos(theta(ia));
		ObservationPoint(2) = z0;
		if (RcsType == 1) //xoz
		{
			ObservationPoint(0) = xory;
			ObservationPoint(1) = 0.0;
			DipoleSFarField(ObservationPoint);
			double ESabsPower2 = 0.0;
			for (int i = 0; i < 3; i++)
			{
				ESabsPower2 = ESabsPower2 + abs(EField(i))*abs(EField(i));
			}
			float Poynting = 0.5*ESabsPower2 / eta;
			radiRCS(ia) = 10 * log10(4 * pi*FarR*FarR*Poynting / P_rad);

		}
		else if (RcsType == 2) //yoz
		{
			ObservationPoint(0) = 0.0;
			ObservationPoint(1) = xory;
			DipoleSFarField(ObservationPoint);
			double ESabsPower2 = 0.0;
			for (int i = 0; i < 3; i++)
			{
				ESabsPower2 = ESabsPower2 + abs(EField(i))*abs(EField(i));
			}
			float Poynting = 0.5*ESabsPower2 / eta;
			radiRCS(ia) = 10 * log10(4 * pi*FarR*FarR*Poynting / P_rad);
		}
		else if (RcsType == 3) //xoy 
		{
			for (int ik = 0; ik < phi.n_elem; ik++)
			{
				ObservationPoint(0) = xory * cos(phi(ik));
				ObservationPoint(1) = xory * sin(phi(ik));
				DipoleSFarField(ObservationPoint);
				double ESabsPower2 = 0.0;
				for (int i = 0; i < 3; i++)
				{
					ESabsPower2 = ESabsPower2 + abs(EField(i))*abs(EField(i));
				}
				float Poynting = 0.5*ESabsPower2 / eta;
				radiRCS(ik) = 10 * log10(4 * pi*FarR*FarR*Poynting / P_rad);
			}
		}
	}
	RadiRCS = radiRCS;
}
