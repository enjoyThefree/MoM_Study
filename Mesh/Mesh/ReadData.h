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
	vector<double> coorData;      // 节点编号对应的坐标
	vector<int> verticesdata;     //组成三角形顶点编号
	Data(int i, vector<double> c) :index(i), coorData(c) {}
	Data(int i, vector<int> v) :index(i), verticesdata(v) {}
};

struct Element {  
	int index;     // element 编号
	int metal_layer_index;
	int dielectric_layer_index;
	//vector<double> pdata;
	vector<int> tdata;     //组成三角形节点编号
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
	int TrianglesNum = 1;    // 三角形数量
	int metal_layer_index;
	int dielectric_layer_index;

	int VerticesNum = 0;
	int ElementNum  = 0;
	bool VerticeCoor_flag;   //是否读取顶点坐标
	bool ElementVerticeNum_flag;       //是否读取 element 顶点编号
	int VerticeNumStart = 0;
	//int ElementStart = 0;
	int VerticeNumbering = 0;
	int ElementNumbering = 0;
	int verticeCount;
	int elementCount;

	string nasFile = "Layer.nas";
	string listFile = "mesh.list";

private:
	//mat t_data;   //element Triangle节点编号信息
	//mat p_data;            // Triangle节点编号的坐标信息

	vector<Data> Vertices_data;   //element Triangle节点编号信息   
	vector<Data> Coordinate_data;   //element Triangle节点编号的坐标信息

	//vector<double> p_data;

	vector<Element> element;

	vector<Layer> Layer_dielectric;   // 介质层信息
	vector<Layer> Layer_metal;        // 金属层信息

};

//class ReadLayerData
//{
//public:
//	void loadLayerinfo();
//	void show_layerInfo();
//private:
//	vector<LayerData> Layer_data;
//};