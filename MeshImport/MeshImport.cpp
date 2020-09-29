#include "MeshImport.h"
#include "../Surface/method.h"
#include <cmath>

namespace Farady
{

std::vector<Element> MeshImport::MeshImportUnit(vector<MeshInfo> &MeshData, vector<FLayer> &metal_layer, 
												vector<FLayer> &dielectric_layer, vector<FPort> &vPort)
{
	for(auto &meshInfo : MeshData){
		MeshInit(meshInfo.filepath);
		MeshRead(metal_layer, dielectric_layer, meshInfo, vPort);

		mergeMesh();

	}

	return Totalelement_list;
}	

void MeshImport::MeshInit(string &filepath)
{
	int pos = filepath.find_last_of('/');
	this->filename = filepath.substr(pos + 1);
	this->filepath = filepath;

    //Determine the file type
    char dot = '.';
	int start = filename.find(dot);	
	filetype = filename.substr(start + 1);

    ifstream infile(filepath.c_str());
	// if (!infile) {
	//     cout << "open " << filename << " file failure !" << endl;
	// 	 exit(0);
    // }

	gridx.clear();
	gridy.clear();
	gridz.clear();
	elementData.clear();
	element_list.clear();

    if(filetype == "nas"){
    	string line;

    	while (infile) {
			getline(infile, line);
			if (infile.fail()) break;

			stringstream sin(line);
			string type;
			sin >> type;
			if (type == "GRID*" || type == "grid*")
				verticesNum++;
			else if (type == "CTRIA3")
				elementNum++;				
    	}
    	gridx.resize(verticesNum + 1);
    	gridy.resize(verticesNum + 1);
    	gridz.resize(verticesNum + 1);
  
      	elementData.resize(elementNum+1);
		
   }
	else if(filetype == "list"){
		string line ,ignore;
		for(int i=0;i<2;i++) getline(infile,line);

    	while (infile) {
			char ch;
			int verticeNum_One;
			int elementNum_One;
			getline(infile, line);
			if (infile.fail()) break;

			stringstream sin(line);
			ch = sin.peek();
			if(ch == '{'){
				getline(infile,line);
				stringstream sin(line);
				sin >> verticeNum_One >>elementNum_One>>ignore;
				verticesNum += verticeNum_One;
				elementNum += elementNum_One;
			}	
		}

		gridx.resize(verticesNum + 1);
    	gridy.resize(verticesNum + 1);
    	gridz.resize(verticesNum + 1);
      	elementData.resize(elementNum+1);
	}

	infile.close();

}

void MeshImport::MeshRead(vector<FLayer> &metal_layer, vector<FLayer> &dielectric_layer, MeshInfo &meshData, vector<FPort> &vPort)
{
	GetCoorOffset(metal_layer, meshData.OnLayer_idx, meshData.x_offset, meshData.y_offset);

	// /*********************  read and save mesh file  ************/

	ifstream infile(filepath.c_str());
	string line ,ignore;
	if(filetype == "nas"){
		while (infile) {
			string type;
			int index1, index2;
			real x, y, z;
			
			getline(infile, line);
			if (infile.fail()) break;

			stringstream sin(line);
			sin >> type;
			if (type == "GRID*" || type == "*") {
				GridLine(line, x, y, z, meshData);
			}
			else if (type == "CTRIA3") {
				CTRIA3Line(line);
			}
			else if (type == "$")
				continue;
			else if (type == "ENDDATA")
				Msg::Debug("EOF!");
			else {
				Msg::Warning("Unknow Input type:");
				Msg::Warning("%s", line.c_str());
			}
		}

	}
	else if(filetype=="list")
	{
		for (int i=0;i<2;i++) getline(infile,line);

		while(infile)
		{
			char ch;
			int verticesNum_One;
			int elementNum_One;

			getline(infile,line);
			stringstream sin(line);

			ch = sin.peek();
			if(ch == '{'){
				VerticeNumStart +=VerticesNum_OneLayer;
				getline(infile,line);
				DataInit(line);
				continue;
			}							
			else if(VerticeCoor_flag == true)
			{
				ReadVerticeCoor(line, meshData);
			}
			else if(ElementVerticeNum_flag == true)
			{
				ReadElementVertice(line);
			}

		}

	}

	infile.close();


	/***************** create new element_list   ******************/
	for (int k = 0; k < elementNum; ++k){
		// element index 从 0 开始
		element_list.push_back(GenerateElement(k, metal_layer, dielectric_layer));
	}

	/*********************  Set Port ************/
	vector<PortElemIdx> vPortElemIdx = meshData.PortElem;
	for(auto &PortElemIdx : vPortElemIdx){
		string PortName = PortElemIdx.PortIdx;
		for(auto &PelemIdx : PortElemIdx.vElemIdx){
			for(auto &elem : element_list){
				if(PelemIdx == elem.GetIndex()){
					// vector<FPoint> rect;
					std::vector<Farady::FPoint> rect = elem.GetVertex();  
					FPort port(rect, portIndex, PortName, true, FDefStd::PortType::PORT_ELEMENT);  //FDefStd::PortType::PORT_ELEMENT
					vPort.push_back(port);
					elem.SetPort(&(vPort.back()));

					++portIndex;
				}
			}
		}
	}

	// CreateElemIdxOnPort_fix(vport);
	// SetPort_fix();


}

void MeshImport::mergeMesh()
{
	for(auto &elem : element_list){
		Totalelement_list.emplace_back(elem);
	}
}

void MeshImport::CreateElemIdxOnPort_fix(vector<FPort> &vport)
{
	ElemIdxOnPort_fix.clear();
	vPort_pointer.clear();
	// for(auto &port : vport){
	for(int j = 0; j < vport.size();++j){
		vector<FPoint> RectPort;
		RectPort = vport[j].GetPortPnts();
		for(auto &elem : element_list){
			vector<FPoint> Vertex = elem.GetVertex();
		
			if(RectPort.size() == Vertex.size()){
				bool IsPort = false;
				for(int i = 0; i < Vertex.size(); ++i){
					if(RectPort[i] == Vertex[i])
						IsPort = true;
					else{
						IsPort = false;
						break;
					}
				}

				if(IsPort){
					ElemIdxOnPort_fix.push_back(elem.GetIndex());

					FPort* port = nullptr;
					port = &vport[j];
					vPort_pointer.push_back(port);
				}
				
			}
		}
	}
}

void MeshImport::GetCoorOffset(vector<FLayer> &metal_layer, int &layer_idx, real &x_offset, real &y_offset)
{
	for(int i = 0; i < metal_layer.size(); i++)
	{
		if(layer_idx == metal_layer[i].GetLayerIndex()){
			this->z_offset = metal_layer[i].GetLayerZStart();
		}
	}
	this->x_offset = x_offset;
	this->y_offset = y_offset;
}


Element MeshImport::GenerateElement(int index, const vector<FLayer>& metal_layer,
	const vector<FLayer>& dielectric_layer)
{
	int i = index + 1;
	int metal_index = 0;
	int diel_index = 0;
	vector<Farady::FPoint> vertex;
	ElementType mesh_type;
	
	double max,min;
	if(elementData[i].size() == 3)
	{
		//vector<FPoint*> elem_pointer;
		mesh_type = TRI;
		FPoint p0(gridx[elementData[i][0]],gridy[elementData[i][0]],gridz[elementData[i][0]]);
    	FPoint p1(gridx[elementData[i][1]],gridy[elementData[i][1]],gridz[elementData[i][1]]);
   	 	FPoint p2(gridx[elementData[i][2]],gridy[elementData[i][2]],gridz[elementData[i][2]]);
		
		// sure element index on port
		//GetPortOnelemIndex(index);

		// CreateElement_points_pointer(p0, p1, p2);
		
		vertex.resize(3);
		vertex[0] = p0;
		vertex[1] = p1;
		vertex[2] = p2;
		
		double z0 = p0[3];
		double z1 = p1[3];
		double z2 = p2[3];
		max = z0 > z1 ? (z0 > z2 ? z0 : z2) : (z1 > z2 ? z1 : z2);    // No need to compare, can be omitted
	    min = z0 < z1 ? (z0 < z2 ? z0 : z2) : (z1 < z2 ? z1 : z2);
	}
	else if(elementData[i].size() == 4){
		//vector<FPoint*> elem_pointer;
		mesh_type = QUAD;
		FPoint p0(gridx[elementData[i][0]],gridy[elementData[i][0]],gridz[elementData[i][0]]);
    	FPoint p1(gridx[elementData[i][1]],gridy[elementData[i][1]],gridz[elementData[i][1]]);
   	 	FPoint p2(gridx[elementData[i][2]],gridy[elementData[i][2]],gridz[elementData[i][2]]);
		FPoint p3(gridx[elementData[i][3]],gridy[elementData[i][3]],gridz[elementData[i][2]]);

		// sure element index on port
		//GetPortOnelemIndex(index);

		// CreateElement_points_pointer(p0, p1, p2, p3);

		vertex.resize(4);
		vertex[0] = p0;
		vertex[1] = p1;
		vertex[2] = p2;
		vertex[3] = p3;

		double z0 = p0[3];
		double z1 = p1[3];
		double z2 = p2[3];
		double z3 = p3[3];
		max = z0 > z1 ? (z0 > z2 ? (z0 > z3 ? z0 : z3) : (z2 > z3 ? z2 : z3)) : (z1 > z2 ? (z1 > z3 ? z1 : z3) : (z2 > z3 ? z2 : z3));
		min = z0 < z1 ? (z0 < z2 ? (z0 < z3 ? z0 : z3) : (z2 < z3 ? z2 : z3)) : (z1 < z2 ? (z1 < z3 ? z1 : z3) : (z2 < z3 ? z2 : z3));
	}
    
	for (const FLayer& m : metal_layer) {
		if (min >= m.GetLayerZStart() && max <= m.GetLayerZStart() + m.GetLayerThickness()) {
			metal_index = m.GetLayerIndex();
			break;
		}
	}
	for (const FLayer& d : dielectric_layer) {
		if (min >= d.GetLayerZStart() && max <= d.GetLayerZStart() + d.GetLayerThickness()) {
			diel_index = d.GetLayerIndex();
			break;
		}
	}

	return Element(index, vertex, metal_index, diel_index, 0,FDefStd::TYPENULL,nullptr, false, mesh_type);
	//return Element(index, vertex, metal_index, diel_index, nullptr);  
}

void MeshImport::SetPort_fix()
{
	if(ElemIdxOnPort_fix.size() == 0){
		Msg::Debug("Warning: The imported object has no port set!"); 
	}
	else{
		for(int i = 0; i < ElemIdxOnPort_fix.size(); ++i){
			for(int j = 0; j < element_list.size(); ++j){
				if(ElemIdxOnPort_fix[i] == j){
					element_list[j].SetPort(vPort_pointer[i]);
				}
			}
		}
	}
	
}

void MeshImport::GridLine(string &line, real &x, real &y, real &z, MeshInfo &meshInfo)
{
	double gain = meshInfo.Gain;
	vector<double> boxMin = meshInfo.BoxMin;

	stringstream sin(line);
	string type;
	int index1, index2;
	//real x, y, z;

	sin >> type;
	if (type == "GRID*") {
		sin >> index1 >> x >> y >> index2;
		assert(index1 == index2);
		x = gain * (x - boxMin.at(0)) + x_offset;
		y = gain * (y - boxMin.at(1)) + y_offset;
		gridx[index1] = x;
		gridy[index1] = y;
	}
	else if (type == "*") {
		sin >> index1 >> z;
		z = gain * z + z_offset;
		// z += z_offset;
		gridz[index1] = z;

	}
}

void MeshImport::GridLineGui(string &line, real &x, real &y, real &z)
{
	stringstream sin(line);
	string type;
	int index1, index2;
	//real x, y, z;

	sin >> type;
	if (type == "GRID*") {
		sin >> index1 >> x >> y >> index2;
		assert(index1 == index2);
		gridx[index1] = x;
		gridy[index1] = y;
	}
	else if (type == "*") {
		sin >> index1 >> z;
		gridz[index1] = z;
		// FPoint p(x,y,z);
		// element_points.push_back(p);

		// cout << "x y z : " << x << " " << y << " " << z << endl;
	}
}

void MeshImport::CTRIA3Line(string &line)
{
	stringstream sin(line);
	string type, ignore;
	int index, p0, p1, p2;
	vector<int> triangle;
	sin >> type >> index >> ignore >> p0 >> p1 >> p2;
    triangle.push_back(p0);
    triangle.push_back(p1);
    triangle.push_back(p2);
	elementData[index] = triangle;
}


void MeshImport::DataInit(string &line)
{
	int verticesNum_One, elementNum_One;
	string ignore;
	stringstream sin(line);
	sin >> verticesNum_One >> elementNum_One >>ignore;
	VerticesNum_OneLayer = verticesNum_One;
	ElementNum_OneLayer = elementNum_One;
	VerticeCoor_flag = true;
	ElementVerticeNum_flag = true;
	verticeCount = 0;
	elementCount = 0;
}

void MeshImport::ReadVerticeCoor(string &line, MeshInfo &meshInfo)
{
	double gain = meshInfo.Gain;
	vector<double> boxMin = meshInfo.BoxMin;

	stringstream sin(line);
	verticeCount++;
	VerticeNumbering++;
	real x, y, z;
	sin >> x >> y >> z;
	x = gain * (x - boxMin.at(0)) + x_offset;
	y = gain * (y - boxMin.at(1)) + y_offset;
	z = gain * z + z_offset;
	// x += x_offset;
	// y += y_offset;
	// z += z_offset;
	gridx[VerticeNumbering] = x;
	gridy[VerticeNumbering] = y;
	gridz[VerticeNumbering] = z;

	// FPoint p(x,y,z);
	// element_points.push_back(p);

	if(verticeCount >= VerticesNum_OneLayer)
		VerticeCoor_flag = false;
}

void MeshImport::ReadVerticeCoorGui(string &line)
{
	stringstream sin(line);

	verticeCount++;
	VerticeNumbering++;
	real x, y, z;
	sin >> x >> y >> z;
	gridx[VerticeNumbering] = x;
	gridy[VerticeNumbering] = y;
	gridz[VerticeNumbering] = z;

	// FPoint p(x,y,z);
	// element_points.push_back(p);

	if(verticeCount >= VerticesNum_OneLayer)
		VerticeCoor_flag = false;
}

void MeshImport::ReadElementVertice(string &line)
{
	stringstream sin(line);
	string ignore;
	elementCount++;
	ElementNumbering++;
	int number;
	int p0, p1, p2, p3;
	sin >> number;
	if(number == 3){
		vector<int> triangle;
		sin >> p0 >> p1 >> p2 >> ignore >> ignore >> ignore;
		p0 +=VerticeNumStart + 1;
		p1 +=VerticeNumStart + 1;
		p2 +=VerticeNumStart + 1;
		triangle.push_back(p0);
		triangle.push_back(p1);
		triangle.push_back(p2);
		elementData[ElementNumbering] = triangle;
		//elementData.push_back(triangle);
	}
	else if(number == 4){
		vector<int> quad;
		sin>>p0>>p1>>p2>>p3>>ignore>>ignore>>ignore;
		p0 +=VerticeNumStart + 1;
		p1 +=VerticeNumStart + 1;
		p2 +=VerticeNumStart + 1;
		p3 +=VerticeNumStart + 1;
		quad.push_back(p0);
		quad.push_back(p1);
		quad.push_back(p2);
		quad.push_back(p3);
		elementData[ElementNumbering] = quad;
		//elementData.push_back(quad);
	}

	if(elementCount >= ElementNum_OneLayer)
		ElementVerticeNum_flag = false;
}



/*****************************************************************    GUI Use display object function         ***********************************/

string MeshImport::GuiDisplay(string &GuifilePath)
{
	MeshInit(GuifilePath);
	MeshGuiRead();
	SaveGUIFile();

	return GuiFileSaveName;
}


void MeshImport::MeshGuiRead()
{
	
	element_list.clear();

	/*********************  read and save mesh file  *************/
	ifstream infile(filepath.c_str());
	string line ,ignore;
	if(filetype == "nas"){
		while (infile) {
			string type;
			int index1, index2;
			real x, y, z;
			
			getline(infile, line);
			if (infile.fail()) break;

			stringstream sin(line);
			sin >> type;
			if (type == "GRID*" || type == "*") {
				GridLineGui(line, x, y, z);
			}
			else if (type == "CTRIA3") {
				CTRIA3Line(line);
			}
			// else if(type == "#" || type == "#Z"){
			// 	GetPortFromFile(line);
			// }
			else if (type == "$")
				continue;
			else if (type == "ENDDATA")
				Msg::Debug("EOF!");
			else {
				Msg::Warning("Unknow Input type:");
				Msg::Warning("%s", line.c_str());
			}
		}

	}
	else if(filetype=="list")
	{
		for (int i=0;i<2;i++) getline(infile,line);

		while(infile)
		{
			char ch;
			int verticesNum_One;
			int elementNum_One;

			getline(infile,line);
			stringstream sin(line);

			ch = sin.peek();
			if(ch == '{'){
				VerticeNumStart +=VerticesNum_OneLayer;
				getline(infile,line);
				DataInit(line);
				continue;
			}							
			else if(VerticeCoor_flag == true)
			{
				ReadVerticeCoorGui(line);
			}
			else if(ElementVerticeNum_flag == true)
			{
				ReadElementVertice(line);
			}

		}

	}

	/***************** create new element_list   ******************/
	for (int k = 0; k < elementNum; ++k){
		// element index 从 0 开始
		element_list.emplace_back(GenerateGuiElement(k));
	}


	infile.close();
}

Element MeshImport::GenerateGuiElement(int index)
{
	int i = index + 1;
	int metal_index = 0;
	int diel_index = 0;
	vector<Farady::FPoint> vertex;
	ElementType mesh_type;
	
	double max,min;
	if(elementData[i].size() == 3)
	{
		mesh_type = TRI;
		FPoint p0(gridx[elementData[i][0]],gridy[elementData[i][0]],gridz[elementData[i][0]]);
    	FPoint p1(gridx[elementData[i][1]],gridy[elementData[i][1]],gridz[elementData[i][1]]);
   	 	FPoint p2(gridx[elementData[i][2]],gridy[elementData[i][2]],gridz[elementData[i][2]]);
		
		vertex.resize(3);
		vertex[0] = p0;
		vertex[1] = p1;
		vertex[2] = p2;
		
	}
	else if(elementData[i].size() == 4){
		mesh_type = QUAD;
		FPoint p0(gridx[elementData[i][0]],gridy[elementData[i][0]],gridz[elementData[i][0]]);
    	FPoint p1(gridx[elementData[i][1]],gridy[elementData[i][1]],gridz[elementData[i][1]]);
   	 	FPoint p2(gridx[elementData[i][2]],gridy[elementData[i][2]],gridz[elementData[i][2]]);
		FPoint p3(gridx[elementData[i][3]],gridy[elementData[i][3]],gridz[elementData[i][2]]);

		vertex.resize(4);
		vertex[0] = p0;
		vertex[1] = p1;
		vertex[2] = p2;
		vertex[3] = p3;

	}
    
	return Element(index, vertex, metal_index, diel_index, 0,FDefStd::TYPENULL,nullptr, false, mesh_type); 
}

void MeshImport::SaveGUIFile()
{
	int pos = filepath.find_last_of('/');
	GuiFileSavePath = filepath.substr(0, pos + 1);
	GuiFileSaveName = GuiFileSavePath + "importNas.mesh";
	char *GuiFileName = (char *)GuiFileSaveName.c_str();

	FILE *fout;
	real Max_x = element_list[0].GetVertex().at(0).XYZ(0);
	real Max_y = element_list[0].GetVertex().at(0).XYZ(1);
	real Max_z = element_list[0].GetVertex().at(0).XYZ(2);
	real Min_x = Max_x;
	real Min_y = Max_y;
	real Min_z = Max_z;
	char outpolyfilename[1000];
	sprintf(outpolyfilename, "%s", GuiFileName);
	fout = fopen(outpolyfilename, "w");
	if (fout == nullptr) {
		Msg::Warning("Can't write mesh_3D output file.");
		return;
	}
	int emi = 0;
	fprintf(fout, "{\n");
	fprintf(fout, "\"vertices\":[");
	for (auto elem : element_list) {
		bool fourelem = false;
		FPoint e4;
		if(elem.GetVertex().size() == 4)
		{
			fourelem = true;
			e4.SetXYZ(elem.GetVertex().at(3).XYZ(0), elem.GetVertex().at(3).XYZ(1), elem.GetVertex().at(3).XYZ(2));
		}
		
		FPoint e1(elem.GetVertex().at(0).XYZ(0), elem.GetVertex().at(0).XYZ(1), elem.GetVertex().at(0).XYZ(2));
		FPoint e2(elem.GetVertex().at(1).XYZ(0), elem.GetVertex().at(1).XYZ(1), elem.GetVertex().at(1).XYZ(2));
		FPoint e3(elem.GetVertex().at(2).XYZ(0), elem.GetVertex().at(2).XYZ(1), elem.GetVertex().at(2).XYZ(2));
		real Ele_Max_X = max(e1.XYZ(0), e2.XYZ(0));
		Ele_Max_X = max(Ele_Max_X, e3.XYZ(0));
		real Ele_Min_X = min(e1.XYZ(0), e2.XYZ(0));
		Ele_Min_X = min(Ele_Min_X, e3.XYZ(0));

		real Ele_Max_Y = max(e1.XYZ(1), e2.XYZ(1));
		Ele_Max_Y = max(Ele_Max_Y, e3.XYZ(1));
		real Ele_Min_Y = min(e1.XYZ(1), e2.XYZ(1));
		Ele_Min_Y = min(Ele_Max_Y, e3.XYZ(1));

		real Ele_Max_Z = max(e1.XYZ(2), e2.XYZ(2));
		Ele_Max_Z = max(Ele_Max_Z, e3.XYZ(2));
		real Ele_Min_Z = min(e1.XYZ(2), e2.XYZ(2));
		Ele_Min_Z = min(Ele_Min_Z, e3.XYZ(2));
		if(fourelem)
		{
			Ele_Max_X = max(Ele_Max_X, e4.XYZ(0));
			Ele_Min_X = min(Ele_Min_X, e4.XYZ(0));
			Ele_Max_Y = max(Ele_Max_Y, e4.XYZ(1));
			Ele_Min_Y = min(Ele_Max_Y, e4.XYZ(1));
			Ele_Max_Z = max(Ele_Max_Z, e4.XYZ(2));
			Ele_Min_Z = min(Ele_Min_Z, e4.XYZ(2));
		}
		if (Max_x <= Ele_Max_X)
			Max_x = Ele_Max_X;
		if (Max_y <= Ele_Max_Y)
			Max_y = Ele_Max_Y;
		if (Max_z <= Ele_Max_Z)
			Max_z = Ele_Max_Z;
		if (Min_x > Ele_Min_X)
			Min_x = Ele_Min_X;
		if (Min_y > Ele_Min_Y)
			Min_y = Ele_Min_Y;
		if (Min_z > Ele_Min_Z)
			Min_z = Ele_Min_Z;
		// cout<<"e1 = "<<e1<<endl;
		// cout<<"e2 = "<<e2<<endl;
		// cout<<"e3 = "<<e3<<endl;
		bool clockwise = true;
		if (e1.XYZ(0) == e2.XYZ(0) && e2.XYZ(0) == e3.XYZ(0)) {
			if (!Method::Ifclockwise(e1.XYZ(2), e1.XYZ(1), e2.XYZ(2), e2.XYZ(1), e3.XYZ(2), e3.XYZ(1)))
				clockwise = false;
		} else if (e1.XYZ(1) == e2.XYZ(1) && e2.XYZ(1) == e3.XYZ(1)) {
			if (!Method::Ifclockwise(e1.XYZ(0), e1.XYZ(2), e2.XYZ(0), e2.XYZ(2), e3.XYZ(0), e3.XYZ(2)))
				clockwise = false;
		} else if (e1.XYZ(2) == e2.XYZ(2) && e2.XYZ(2) == e3.XYZ(2)) {
			// cout<<"3"<<endl;
			if (!Method::Ifclockwise(e1.XYZ(0), e1.XYZ(1), e2.XYZ(0), e2.XYZ(1), e3.XYZ(0), e3.XYZ(1)))
				clockwise = false;
		} else {
			clockwise = false;
			// cout<<endl;
			// cout<<"e1 = "<<e1<<endl;
			// cout<<"e2 = "<<e2<<endl;
			// cout<<"e3 = "<<e3<<endl;
			// // cout<<"4"<<endl;
			// // cout<<"e1.XYZ(2) = "<<e1.XYZ(2)<<"e2.XYZ(2) = "<<e2.XYZ(2)<<"e3.XYZ(2) = "<<e3.XYZ(2)<<endl;
			// assert(0);
		}
		if (emi == element_list.size() - 1) {
				if(fourelem)
				fprintf(fout, "[%f,%f,%f],[%f,%f,%f],[%f,%f,%f],[%f,%f,%f]", e1.XYZ(0), e1.XYZ(1), e1.XYZ(2), e2.XYZ(0), e2.XYZ(1),
						e2.XYZ(2), e3.XYZ(0), e3.XYZ(1), e3.XYZ(2), e4.XYZ(0), e4.XYZ(1), e4.XYZ(2));				
				else
				fprintf(fout, "[%f,%f,%f],[%f,%f,%f],[%f,%f,%f]", e1.XYZ(0), e1.XYZ(1), e1.XYZ(2), e2.XYZ(0), e2.XYZ(1),
						e2.XYZ(2), e3.XYZ(0), e3.XYZ(1), e3.XYZ(2));
			// OutMeshdata<<"["<<e3.XYZ(0)<<","<<e3.XYZ(1)<<","<<e3.XYZ(2)<<"],"<<"["<<e2.XYZ(0)<<","<<e2.XYZ(1)<<","<<e2.XYZ(2)<<"],"<<"["<<e1.XYZ(0)<<","<<e1.XYZ(1)<<","<<e1.XYZ(2)<<"]";
		} else {
				if (fourelem)
				{
					// cout<<"e4 = "<<e4<<"e3 = "<<e3<<endl;
					fprintf(fout, "[%f,%f,%f],[%f,%f,%f],[%f,%f,%f],[%f,%f,%f],", e1.XYZ(0), e1.XYZ(1), e1.XYZ(2), e2.XYZ(0), e2.XYZ(1),
						e2.XYZ(2), e3.XYZ(0), e3.XYZ(1), e3.XYZ(2), e4.XYZ(0), e4.XYZ(1), e4.XYZ(2));	
				}				
				else
				fprintf(fout, "[%f,%f,%f],[%f,%f,%f],[%f,%f,%f],", e1.XYZ(0), e1.XYZ(1), e1.XYZ(2), e2.XYZ(0),
						e2.XYZ(1), e2.XYZ(2), e3.XYZ(0), e3.XYZ(1), e3.XYZ(2));
			// if (clockwise)
			// {
			// 	if (fourelem)
			// 	{
			// 		// cout<<"e4 = "<<e4<<"e3 = "<<e3<<endl;
			// 		fprintf(fout, "[%f,%f,%f],[%f,%f,%f],[%f,%f,%f],[%f,%f,%f],", e1.XYZ(0), e1.XYZ(1), e1.XYZ(2), e2.XYZ(0), e2.XYZ(1),
			// 			e2.XYZ(2), e3.XYZ(0), e3.XYZ(1), e3.XYZ(2), e4.XYZ(0), e4.XYZ(1), e4.XYZ(2));	
			// 	}				
			// 	else
			// 	fprintf(fout, "[%f,%f,%f],[%f,%f,%f],[%f,%f,%f],", e1.XYZ(0), e1.XYZ(1), e1.XYZ(2), e2.XYZ(0),
			// 			e2.XYZ(1), e2.XYZ(2), e3.XYZ(0), e3.XYZ(1), e3.XYZ(2));
			// }
			// // OutMeshdata<<"["<<e1.XYZ(0)<<","<<e1.XYZ(1)<<","<<e1.XYZ(2)<<"],"<<"["<<e2.XYZ(0)<<","<<e2.XYZ(1)<<","<<e2.XYZ(2)<<"],"<<"["<<e3.XYZ(0)<<","<<e3.XYZ(1)<<","<<e3.XYZ(2)<<"],";
			// else
			// {
			// 	if (fourelem)
			// 	{
			// 		// cout<<"e4 = "<<e4<<"e3 = "<<e3<<endl;
			// 		fprintf(fout, "[%f,%f,%f],[%f,%f,%f],[%f,%f,%f],[%f,%f,%f],", e4.XYZ(0), e4.XYZ(1), e4.XYZ(2),
			// 				e3.XYZ(0), e3.XYZ(1), e3.XYZ(2), e2.XYZ(0), e2.XYZ(1), e2.XYZ(2), e1.XYZ(0), e1.XYZ(1),
			// 				e1.XYZ(2));					
			// 	}
			// 	else
			// 		fprintf(fout, "[%f,%f,%f],[%f,%f,%f],[%f,%f,%f],", e3.XYZ(0), e3.XYZ(1), e3.XYZ(2), e2.XYZ(0),
			// 				e2.XYZ(1), e2.XYZ(2), e1.XYZ(0), e1.XYZ(1), e1.XYZ(2));
			// }
			// OutMeshdata<<"["<<e3.XYZ(0)<<","<<e3.XYZ(1)<<","<<e3.XYZ(2)<<"],"<<"["<<e2.XYZ(0)<<","<<e2.XYZ(1)<<","<<e2.XYZ(2)<<"],"<<"["<<e1.XYZ(0)<<","<<e1.XYZ(1)<<","<<e1.XYZ(2)<<"],";
		}
		// cout<<"["<<e1.XYZ(0)<<","<<elem.XYZ(1)<<","<<elem.XYZ(2)<<"],";
		emi++;
	}
	fprintf(fout, "],\n");
	int eii = 0;
	emi = 0;
	fprintf(fout, "\"faces\":[");
	for (auto elem : element_list) {
		if (emi == element_list.size() - 1) {
			if (elem.GetVertex().size() == 4) {
				fprintf(fout, "[%d,%d,%d,%d,100]", eii, eii + 1, eii + 2, eii + 3);
				eii = eii + 4;
			} else {
				fprintf(fout, "[%d,%d,%d,100]", eii, eii + 1, eii + 2);
				eii = eii + 3;
			}

		}
		// OutMeshdata<<"["<<eii<<","<<eii+1<<","<<eii+2<<",100"<<"]";
		else {
			if (elem.GetVertex().size() == 4) {
				fprintf(fout, "[%d,%d,%d,%d,100],", eii, eii + 1, eii + 2, eii + 3);
				eii = eii + 4;
			} else {
				fprintf(fout, "[%d,%d,%d,100],", eii, eii + 1, eii + 2);
				eii = eii + 3;
			}
		}
		// OutMeshdata<<"["<<eii<<","<<eii+1<<","<<eii+2<<",100"<<"],";
		emi++;
		
	}
	
	real L_x = Max_x - Min_x;
	real W_y = Max_y - Min_y;
	real H_z = Max_z - Min_z;
	real Big_LWH = max(L_x, W_y);
	Big_LWH = max(H_z, Big_LWH);
	fprintf(fout, "],\n");
	// fprintf(fout, "\"bbox\":[%f,%f,%f,%f]\n", Max_x, Max_y, Max_z, Big_LWH);   //out mesh add Min_x,Min_y,Min_z
	fprintf(fout, "\"bbox\":[%f,%f,%f,%f,%f,%f,%f]\n", Max_x, Max_y, Max_z, Min_x, Min_y, Min_z, Big_LWH);
	fprintf(fout, "}\n");
	fclose(fout);
}


/*****************************************************************    python get box min function         ***********************************/
vector<double> MeshImport::PythonGetBoxMin(string &filepath)
{
	double x_min, y_min;

	int pos = filepath.find_last_of('/');
	this->filename = filepath.substr(pos + 1);
	this->filepath = filepath;

    //Determine the file type
    char dot = '.';
	int start = filename.find(dot);	
	filetype = filename.substr(start + 1);

    ifstream infile(filepath.c_str());
	string line ,ignore;
	if(filetype == "nas"){
		while (infile) {
			string type;
			int index1, index2;
			real x, y, z;
			
			getline(infile, line);
			if (infile.fail()) break;

			stringstream sin(line);
			sin >> type;
			if (type == "GRID*") {
				sin >> index1 >> x >> y >> index2;
				assert(index1 == index2);
			
				if(index1 == 1 || (x < x_min && y < y_min)){
					x_min = x;
					y_min = y;
				}
				else{
					continue;
				}
			
			}
			else if(type == "*"){
				continue;
			}
			else if (type == "CTRIA3") {
				continue;
			}
			else if (type == "$")
				continue;
			else if (type == "ENDDATA")
				Msg::Debug("EOF!");
			else {
				Msg::Warning("Unknow Input type:");
				Msg::Warning("%s", line.c_str());
			}
		}

	}
	else if(filetype=="list")
	{

		for (int i=0;i<2;i++) getline(infile,line);

		while(infile)
		{
			char ch;
			int verticesNum_One;
			int elementNum_One;

			getline(infile,line);
			stringstream sin(line);

			ch = sin.peek();
			if(ch == '{'){
				VerticeNumStart +=VerticesNum_OneLayer;
				getline(infile,line);
				DataInit(line);
				continue;
			}							
			else if(VerticeCoor_flag == true)
			{
				verticeCount++;
				VerticeNumbering++;
				real x, y, z;
				sin >> x >> y >> z;

				if(VerticeNumbering == 1 || (x < x_min && y < y_min)){
					x_min = x;
					y_min = y;
				}
				else{
					continue;
				}
		

				if(verticeCount >= VerticesNum_OneLayer)
					VerticeCoor_flag = false;
			}
			else if(ElementVerticeNum_flag == true)
			{
				continue;
			}

		}

	}

	vector<double> BoxMin;
	BoxMin.resize(2);
	BoxMin[0] = x_min;
	BoxMin[1] = y_min;

	return BoxMin;
}

}   // namespace Faraday