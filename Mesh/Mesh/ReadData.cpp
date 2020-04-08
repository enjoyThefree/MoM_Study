/**************************   ReadData.cpp     **************************/
#include"ReadData.h"

bool ReadData::loadData()
{
	Vertices_data.clear();
	Coordinate_data.clear();
	element.clear();
	loadLayerinfo();  //提供 Layer 信息

	//string nasFile = "test.nas";
	ifstream infile(nasFile.c_str());
	if (!infile) {
		cout << "open nas file failure !" << endl;
		exit(0);
	}

	string line;
	for (int i = 0; i < 6; i++) getline(infile, line);

	//get TrianglesNum
	string str;
	int trianglesNum;
	getline(infile, line);
	stringstream sin(line);
	while (sin) {
		sin >> str >> str;
		sin >> trianglesNum;
		TrianglesNum = trianglesNum;
	}

	for(int i=0;i<2;i++)  getline(infile, line);

	while (infile) {
		string ignore;
		int index;
		double x, y, z;
		int p1, p2, p3;
		vector<double> coordata;
		vector<int> verticedata;
		getline(infile, line);

		if (infile.fail()) break;     //防止最后一行数据读两遍

		stringstream sin(line);
		sin >> str;
		if (str == "GRID*" || str == "*") {  //节点编号的坐标
			if (str == "GRID*") {
				sin >> index >> x >> y >> ignore;
				/*T_data.push_back(x);
				T_data.push_back(y);*/
			}
			else if (str == "*") {
				sin >> index >> z;
				coordata.push_back(x);
				coordata.push_back(y);
				coordata.push_back(z);
				Coordinate_data.push_back(Data(index, coordata));
			}
			
		}
		else if (str == "CTRIA3") { //组成三角形节点编号
			sin >> index >> ignore >> p1 >> p2 >> p3;

			//给三个点的z坐标排序，取最大和最小的
			double z1 = Coordinate_data[p1 - 1].coorData[2];
			double z2 = Coordinate_data[p2 - 1].coorData[2];
			double z3 = Coordinate_data[p3 - 1].coorData[2];
			double max = z1 > z2 ? (z1 > z3 ? z1 : z3) : (z2 > z3 ? z2 : z3);
			double min = z1 < z2 ? (z1 < z3 ? z1 : z3) : (z2 < z3 ? z2 : z3);
			for (int i = 0; i < Layer_dielectric.size(); i++) {
				if (min >= Layer_dielectric[i].data.first && max <= (Layer_dielectric[i].data.first + Layer_dielectric[i].data.second)) {
					dielectric_layer_index = Layer_dielectric[i].index;
				}
			}
			for (int i = 0; i < Layer_metal.size(); i++) {
				if (min >= Layer_metal[i].data.first && max <= (Layer_metal[i].data.first + Layer_metal[i].data.second)) {
					metal_layer_index = Layer_metal[i].index;
				}
			}
			verticedata.push_back(p1);
			verticedata.push_back(p2);
			verticedata.push_back(p3);
			Vertices_data.push_back(Data(index, verticedata));
			element.push_back(Element(index, metal_layer_index, dielectric_layer_index, verticedata));
		}
	}

	infile.close();

	return true;
}

void ReadData::loadListData()
{
	ifstream infile(listFile.c_str());
	if (!infile.is_open()) {
		cout << "open list file failure !" << endl;
		exit(0);
	}

	string line;
	string ignore;
	vector<double> coordata;
	vector<int> verticedata;

	for (int i = 0; i < 2; i++) getline(infile, line);
	
	while (infile) {
		int verticesNum;
		int elementNum;
		/*int verticeCount;
		int elementCount;*/
		char ch;
		string str;
		getline(infile, line);
		stringstream sin(line);
		ch = sin.peek();
		if (ch == '{') {
			VerticeNumStart += VerticesNum;
			getline(infile, line);
			stringstream sin(line);
			sin >> verticesNum >> elementNum >> ignore;
			VerticesNum = verticesNum;
			ElementNum = elementNum;
			VerticeCoor_flag = true;
			ElementVerticeNum_flag = true;
			verticeCount = 0;
			elementCount = 0;
			continue;
		}	
		
		if (VerticeCoor_flag == true) {  //顶点坐标
			coordata.clear();
			verticeCount++;
			VerticeNumbering++;
			double x, y, z;
			sin >> x >> y >> z;
			coordata.push_back(x);
			coordata.push_back(y);
			coordata.push_back(z);
			Coordinate_data.push_back(Data(VerticeNumbering,coordata));   //顶点编号，对应坐标
			
			if (verticeCount >= VerticesNum)
				VerticeCoor_flag = false;
		}
		else if(ElementVerticeNum_flag == true)   
		{
			verticedata.clear();
			elementCount++;
			ElementNumbering++;
			int number;
			int p1, p2, p3, p4;
			sin >> number;
			if (number == 3) {
				sin >> p1 >> p2 >> p3 >> ignore >> ignore >> ignore;
				p1 += VerticeNumStart+1;
				p2 += VerticeNumStart+1;
				p3 += VerticeNumStart+1;
				verticedata.push_back(p1);
				verticedata.push_back(p2);
				verticedata.push_back(p3);
				Vertices_data.push_back(Data(ElementNumbering, verticedata));    //element的编号，组成element的顶点编号

				//判断介质层、金属层
				double z1 = Coordinate_data[p1 - 1].coorData[2];
				double z2 = Coordinate_data[p2 - 1].coorData[2];
				double z3 = Coordinate_data[p3 - 1].coorData[2];
				double max = z1 > z2 ? (z1 > z3 ? z1 : z3) : (z2 > z3 ? z2 : z3);
				double min = z1 < z2 ? (z1 < z3 ? z1 : z3) : (z2 < z3 ? z2 : z3);
				for (int i = 0; i < Layer_dielectric.size(); i++) {
					if(min >=Layer_dielectric[i].data.first && max<=(Layer_dielectric[i].data.first+Layer_dielectric[i].data.second))
						dielectric_layer_index = Layer_dielectric[i].index;
				}
				for (int i = 0; i < Layer_metal.size(); i++) {
					if (min >= Layer_metal[i].data.first && max <= (Layer_metal[i].data.first + Layer_metal[i].data.second))
						metal_layer_index = Layer_metal[i].index;
				}
				element.push_back(Element(ElementNumbering, metal_layer_index, dielectric_layer_index, verticedata));
			}
			else if (number == 4)
			{
				sin >> p1 >> p2 >> p3 >> p4 >> ignore >> ignore >> ignore;
				p1 += VerticeNumStart + 1;
				p2 += VerticeNumStart + 1;
				p3 += VerticeNumStart + 1;
				p4 += VerticeNumStart + 1;
				verticedata.push_back(p1);
				verticedata.push_back(p2);
				verticedata.push_back(p3);
				verticedata.push_back(p4);
				Vertices_data.push_back(Data(ElementNumbering, verticedata));    //element的编号，组成element的顶点编号

				double z1 = Coordinate_data[p1 - 1].coorData[2];
				double z2 = Coordinate_data[p2 - 1].coorData[2];
				double z3 = Coordinate_data[p3 - 1].coorData[2];
				double z4 = Coordinate_data[p4 - 1].coorData[2];
				double max = z1 > z2 ? (z1 > z3 ? (z1 > z4 ? z1 : z4) : (z3 > z4 ? z3 : z4)) : (z2 > z3 ? (z2 > z4 ? z2 : z4) : (z3 > z4 ? z3 : z4));
				double min = z1 < z2 ? (z1 < z3 ? (z1 < z4 ? z1 : z4) : (z3 < z4 ? z3 : z4)) : (z2 < z3 ? (z2 < z4 ? z2 : z4) : (z3 < z4 ? z3 : z4));
				for (int i = 0; i < Layer_dielectric.size(); i++) {
					if (min >= Layer_dielectric[i].data.first && max <= (Layer_dielectric[i].data.first + Layer_dielectric[i].data.second))
						dielectric_layer_index = Layer_dielectric[i].index;
				}
				for (int i = 0; i < Layer_metal.size(); i++) {
					if (min >= Layer_metal[i].data.first && max <= (Layer_metal[i].data.first + Layer_metal[i].data.second))
						metal_layer_index = Layer_metal[i].index;
				}
				element.push_back(Element(ElementNumbering, metal_layer_index, dielectric_layer_index, verticedata));
			}

			if (elementCount > ElementNum)
				ElementVerticeNum_flag = false;
		}

	}
	infile.close();
}

void ReadData::loadLayerinfo() {  //用于测试用，提供 Layer信息。
	Layer_dielectric.clear();
	Layer_metal.clear();

	ifstream infile("Layer.txt");
	if (!infile) {
		cout << "open Layer file failure !" << endl;
		exit(0);
	}

	string line;
	/*getline(infile, line);*/

	while (infile) {
		string str;
		int index;
		//char ch;
		double start, thickness;
		pair<double,double> data;
		getline(infile, line);

	    if (infile.fail()) break;     //防止最后一行数据读两遍

		stringstream sin(line);
		//ch = sin.peek();
		sin >> str;
		if (str == "dielectric") {
			sin >> index >> start >> thickness;
			data.first = start;
			data.second = thickness;
			Layer_dielectric.push_back(Layer(index, data));
		}
		else if (str == "metal") {
			sin >> index >> start >> thickness;
			data.first = start;
			data.second = thickness;
			Layer_metal.push_back(Layer(index, data));
		}
		
	}

	infile.close();
}

void ReadData::show()
{
	ofstream outf1("Coordinate_data.txt");
	ofstream outf2("Vertices_data.txt");

	if (!outf1.is_open()) {
		cout << "open p_data file failure !" << endl;
		exit(0);
	}

	if (!outf2.is_open()) {
		cout << "open t_data file failure !" << endl;
		exit(0);
	}

	for (int i = 0; i < Coordinate_data.size(); i++) {
		outf1 << Coordinate_data[i].coorData[0] << " " << Coordinate_data[i].coorData[1] << " " << Coordinate_data[i].coorData[2] << endl;
	}

	for (int i = 0; i < Vertices_data.size(); i++) {
		if(Vertices_data[i].verticesdata.size() == 3)
			outf2 << Vertices_data[i].verticesdata[0] << " " << Vertices_data[i].verticesdata[1] << " " << Vertices_data[i].verticesdata[2] << endl;
		else if(Vertices_data[i].verticesdata.size() == 4 )
			outf2 << Vertices_data[i].verticesdata[0] << " " << Vertices_data[i].verticesdata[1] << " " << Vertices_data[i].verticesdata[2] << " " << Vertices_data[i].verticesdata[3] << endl;
	}

	cout << "保存完毕，退出中..." << endl;
	outf1.close();
	outf2.close();

}

void ReadData::show_layerInfo() {
	for (int i = 0; i < Layer_dielectric.size(); i++) {
		cout << "编号：" << Layer_dielectric[i].index << "\t" << "起点：" << Layer_dielectric[i].data.first << "\t" << "厚度：" << Layer_dielectric[i].data.second << endl;
	}
	for (int i = 0; i < Layer_metal.size(); i++) {
		cout << "编号：" << Layer_metal[i].index << "\t" << "起点：" << Layer_metal[i].data.first << "\t" << "厚度：" << Layer_metal[i].data.second << endl;
	}
}