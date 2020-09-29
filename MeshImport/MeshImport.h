#ifndef MeshImport_H
#define MeshImport_H

#include "../../Geomtry/GeomStdHeads.h"
#include "../Global/Element.h"
#include "../../Geomtry/FGeomBaseStd/FPort.h"
#include "../../Geomtry/FBaseStd/FLayer.h"
#include "../../Common/EMUtils.h"
#include<vector>

using namespace std;


namespace Farady
{

class PortElemIdx
{
public:
	PortElemIdx(){}
	~PortElemIdx(){}

	PortElemIdx(string s, vector<int> elemIdx){
		PortIdx = s;
		vElemIdx = elemIdx;
	}

	bool operator==(const PortElemIdx& PEIdx){return false;}
	PortElemIdx(const PortElemIdx& PEIdx)
	{
	    PortIdx = PEIdx.PortIdx;
		vElemIdx = PEIdx.vElemIdx;
	}
	const PortElemIdx& operator=(const PortElemIdx& PEIdx)
	{
	    if(this != &PEIdx){
	        PortIdx = PEIdx.PortIdx;
			vElemIdx = PEIdx.vElemIdx;
	    }
        return *this;
	}

public:
	string PortIdx;
	vector<int> vElemIdx;

};

class MeshInfo
{
public:
	MeshInfo(){}

	MeshInfo(int &OnLayer_idx, double &x_offset, double &y_offset, string &filepath);

	~MeshInfo(){}

	bool operator==(const MeshInfo& mi){return false;}
	MeshInfo(const MeshInfo& mi)
	{
	    OnLayer_idx = mi.OnLayer_idx;
	    x_offset = mi.x_offset;
	    y_offset = mi.y_offset;
	    filepath = mi.filepath;
		Gain = mi.Gain;
		BoxMin = mi.BoxMin;
		PortElem = mi.PortElem;
	}
	const MeshInfo& operator=(const MeshInfo& mi)
	{
	    if(this != &mi){
	        OnLayer_idx = mi.OnLayer_idx;
            x_offset = mi.x_offset;
            y_offset = mi.y_offset;
            filepath = mi.filepath;
			Gain = mi.Gain;
			BoxMin = mi.BoxMin;
			PortElem = mi.PortElem;
	    }
        return *this;
	}

public:
	int OnLayer_idx;
	double x_offset;
	double y_offset;
	string filepath;

	double Gain = 1;               // Import mesh file magnification
	vector<double> BoxMin;	       // Boundarying box min Value.   0:x_min ; 1:y_min
	vector<PortElemIdx> PortElem;  

};

class MeshImport
{
public:
	MeshImport(){}

	~MeshImport(){
		vPort_pointer.clear();
		gridx.clear();
		gridy.clear();
		gridz.clear();
		elementData.clear();
		element_list.clear();
		ElemIdxOnPort_fix.clear();
		Totalelement_list.clear();
	}
	
	std::vector<Element> MeshImportUnit(vector<MeshInfo> &MeshData, vector<FLayer> &metal_layer, vector<FLayer> &dielectric_layer, vector<FPort> &vPort);

	vector<Element> GetTotalElement() { return this-> Totalelement_list; }

	// GUI Display object
	string GuiDisplay(string &GuifilePath);

	// python get box min
	vector<double> PythonGetBoxMin(string &filepath);
	
private:
	std::vector<Element> Totalelement_list;
	std::vector<Element> element_list;
	string filepath;
	real x_offset = 0, y_offset = 0, z_offset = 0;

	string filename;
    string filetype;
	string GuiFileSavePath;
	string GuiFileSaveName;

	int verticesNum =0;
	int elementNum = 0;

	int VerticesNum_OneLayer = 0;
	int ElementNum_OneLayer = 0;
	bool VerticeCoor_flag;         //Whether to read vertex coordinates
	bool ElementVerticeNum_flag;   //Whether to read the element vertex number
	int VerticeNumStart = 0;
	int VerticeNumbering = 0;
	int ElementNumbering = 0;
	int verticeCount;
	int elementCount;

	int portIndex = 0;

private:

	vector<double> gridx, gridy, gridz;
    vector<vector<int>> elementData;

	vector<int> ElemIdxOnPort_fix;
	vector<FPort*> vPort_pointer;

private:
	void MeshInit(string &filepath);
	
	void MeshRead(vector<FLayer> &metal_layer, vector<FLayer> &dielectric_layer, 
					MeshInfo &MeshData, vector<FPort> &vPort);

	void mergeMesh();

	void CreateElemIdxOnPort_fix(vector<FPort> &vport);

	Element GenerateElement(int index, const vector<FLayer>& metal_layer,
		const vector<FLayer>& dielectric_layer);
	
	void GetCoorOffset(vector<FLayer> &metal_layer, int &layer_idx, real &x_offset, real &y_offset);


	void SetPort_fix();

	//Read Mesh NasFile
	void GridLine(string &line, real &x, real &y, real &z, MeshInfo &meshInfo);
	void GridLineGui(string &line, real &x, real &y, real &z);
	void CTRIA3Line(string &line);

	// Read Mesh ListFile
	void DataInit(string &line);
	void ReadVerticeCoor(string &line, MeshInfo &meshInfo);
	void ReadVerticeCoorGui(string &line);		// Gui use
	void ReadElementVertice(string &line);

	// GUI Use display object
	void MeshGuiRead();
	Element GenerateGuiElement(int index);
	void SaveGUIFile();

};

} //namespace Faraday

#endif // MeshImport_H