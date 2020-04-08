#include<iostream>
#include<fstream>
#include<armadillo>
#include"ReadData.h"

using namespace std;
using namespace arma;

int main()
{
	ReadData read;
	read.loadListData();
	/*read.loadData();
	read.loadLayerinfo();
	read.show_layerInfo();*/
	read.show();

	cout << "Hello World !" << endl;
}

