//#pragma once
//#include"InputData.h"
//
//namespace mom {
//
//	bool InputStream(ifstream& in, stringstream& ss) {
//		string str;
//		ss.clear();
//		getline(in, str);
//		ss << str.substr(0, str.find("#"));
//		return true;
//	}
//
//	bool InputData::readData(mat &basis_list, cube &element_list) {
//		string ignore;
//		stringstream ss;
//		/*读取rwg:个数以及相邻的两个element*/
//		int num = 0, basis_p = 0, basis_n = 0;
//		ifstream basis_in("rwg.txt",ios::in);
//		if (!basis_in.is_open())
//			cout << "open rwg.txt file failure";
//
//		//basis_in >> num;
//		InputStream(basis_in,ss);
//		ss >> num;
//		basis_list.set_size(num, 2);
//		for (int k = 0; k < num; k++) {
//			InputStream(basis_in, ss);
//			ss >> basis_p >> basis_n >> ignore >> ignore;
//			basis_list(k, 0) = basis_p;
//			basis_list(k, 1) = basis_n;
//		}
//		//int k;
//		//while (num--) {
//		//	//basis_in >> basis_p >> basis_n >> ignore >> ignore;
//		//	InputStream(basis_in, ss);
//		//	ss >> basis_p >> basis_n >> ignore >> ignore;
//		//	basis_list(k, 0) = basis_p;
//		//	basis_list(k, 1) = basis_n;
//		//	k++;
//		//}
//		basis_in.close();
//
//		/*读取坐标,element编号*/
//		ifstream element_in("tri.txt", ios::in);
//		double x, y, z;
//		//mat vertex_list(1, 3);
//		if (!element_in.is_open())
//			cout << "open tri.txt file failure";
//		InputStream(element_in, ss);
//		ss >> num;
//		element_list.set_size(3, 3,num);  //第一层为第一个element
//		for (int i = 0; i < num; i++) {
//			mat tem(3, 3);  //第一行：p1点 ...  类推
//			for (int j = 0; j < 3; j++) {
//				InputStream(element_in, ss);
//				ss >> x >> y >> z;
//				tem(j, 0) = x;
//				tem(j, 1) = y;
//				tem(j, 2) = z;
//			}
//			element_list.slice(i) = tem;
//		}
//		element_in.close();
//
//		return true;
//	}
//
//}