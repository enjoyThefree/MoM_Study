/*********************    ReadData.h   **********************/
#pragma once
#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<sstream>
#include<armadillo>

using namespace std;
using namespace arma;

//struct Data {
//	int index;
//	vector<double> data;
//	Data(int i, vector<double> t) :index(i), data(t) {}
//};

struct Data {
	int index;     
	int metal_layer_index;
	int dielectric_layer_index;
	vector<double> coorData;      // �ڵ��Ŷ�Ӧ������
	vector<int> verticesdata;     //��������ζ�����
	Data(int i, vector<double> c) :index(i), coorData(c) {}
	Data(int i, vector<int> v) :index(i), verticesdata(v) {}
};

struct Element {  
	int index;     // element ���
	int metal_layer_index;
	int dielectric_layer_index;
	//vector<double> pdata;
	vector<int> tdata;     //��������νڵ���
	//Data(int i, vector<double> t) :index(i), pdata(t) {}
	//Data(int i, vector<int> t) :index(i), tdata(t) {}
	//Data(int i,int metal_index,int diel_index):index(i),metal_layer_index(metal_index),dielectric_layer_index(diel_index) {}
	Element(int i,int metal_index,int diel_index,vector<int> t):index(i),metal_layer_index(metal_index),dielectric_layer_index(diel_index),tdata(t) {}
};

struct Layer { 
	int index;
	pair<double,double> data;   // first is start and second is thickness
	Layer(int i, pair<double, double> d):index(i),data(d) {}
};

struct ListData {
	
};

class ReadData 
{
public:
	bool loadData();
	void loadLayerinfo();
	void loadListData();
	void show();
	void show_layerInfo();
	//void show_metal_diel_index();

private:
	int TrianglesNum = 1;    // ����������
	int metal_layer_index;
	int dielectric_layer_index;

	int VerticesNum = 0;
	int ElementNum  = 0;
	bool VerticeCoor_flag;   //�Ƿ��ȡ��������
	bool ElementVerticeNum_flag;       //�Ƿ��ȡ element ������
	int VerticeNumStart = 0;
	//int ElementStart = 0;
	int VerticeNumbering = 0;
	int ElementNumbering = 0;
	int verticeCount;
	int elementCount;

	string nasFile = "Layer.nas";
	string listFile = "mesh.list";

private:
	//mat t_data;   //element Triangle�ڵ�����Ϣ
	//mat p_data;            // Triangle�ڵ��ŵ�������Ϣ

	vector<Data> Vertices_data;   //element Triangle�ڵ�����Ϣ   
	vector<Data> Coordinate_data;   //element Triangle�ڵ��ŵ�������Ϣ

	//vector<double> p_data;

	vector<Element> element;

	vector<Layer> Layer_dielectric;   // ���ʲ���Ϣ
	vector<Layer> Layer_metal;        // ��������Ϣ

};

//class ReadLayerData
//{
//public:
//	void loadLayerinfo();
//	void show_layerInfo();
//private:
//	vector<LayerData> Layer_data;
//};