#pragma once
#include<iostream>
#include<vector>
#include<armadillo>
#include <sstream>

using namespace std;
using namespace arma;

class FPort;
class FPoint;

enum MeshType {TRIANGLE=0,TETRAHEDRON};

class Element
{
private:
	int index;  //index of element

	vector<FPoint> element_vertex;   //element vertex 的坐标

	vector<vec3> element_vertex_vec;

	int metal_layer_index;
	int dielectric_layer_index;
	MeshType mesh_type;

	FPort* port;

	bool is_TDS_interface;    //

	bool is_TDS;

	vector<double>gauss_weight;
	vector<double> diff_gauss_weight;

	vector<vec3> gauss_points_vec;
	vector<vec3> diff_gauss_points_vec;

	static bool init_flag;
	static double sym_gauss_weight_tet[9][100];
	static double sym_gauss_coord_tet[9][4][100];
	static int num_sym_gauss_nodes_tet[9];
	static void GenerateWeightaAndCoordinatesForSymGaussianQuadrature();

public:
	Element(MeshType mesh_type = TRIANGLE);
	Element(int index, const vector<FPoint>& element_vertex, int metal_layer_index,
		int dielectric_layer_index, FPoint* ppoint = nullptr,
		bool is_TDS_interface, MeshType mesh_type = TRIANGLE);

	void SetTDSInterface(bool is_TDS_interface) { this->is_TDS_interface(); }
	bool IsTDSInterface() const { return this->is_TDS_interface; }

	void SetTDS(bool is_TDS) { this->is_TDS = is_TDS; }
	bool IsTDS()const { return this->is_TDS; }

	bool isTDSElement() const { return this->IsTDS() || IsTDSInterface(); }
	~Element();

	void GenerateGauss();
	void Initialize(MeshType mesh_type = TRIANGLE);

	vector<bool> visited;

	//邻接的三角形列表
	vector<pair<int, int>> adjacent_info_ix;

	vector<int> shared_metal_layer_idxes;

	//Get the values of the member variables
	int GetIndex() const;
	void SetIndex(int SetIndex) { this->index = SetIndex; }
	const vector<FPoint>& GetVertex() const;
	vector<FPoint>& GetVertex();

	const vector<vec3>& GetVertexVec() const { return this->element_vertex_vec; }

	FPoint* GetPort() const;

	//return true if this element is on port
	FPoint GetMiddlePoint() const;
	void SetPort(FPort* port);
	int GetMetalLayerIndex() const;
	void SetDielectricLayerIndex(int die_layer_idx) { this->dielectric_layer_index = die_layer_idx; }

	const vector<double> &GetGauseWeight() const { return gauss_weight; }
	const vector<double> &GetDiffGauseWeight() const { return diff_gauss_weight; }
	const vector<vec3>& GetGaussPointVec() const { return gauss_points_vec; }
	const vector<vec3>& GetDiffGaussPointVec() const { return diff_gauss_points_vec; }
	const vector<double>* GetGaussWeightP() const { return &gauss_weight; }
	const vector<double>* GetDiffGaussWeightP() const { return &diff_gauss_weight; }
	const vector<vec3>* GetGaussPointVecP() const { return &gauss_points_vec; }
	const vector<vec3>* GetDiffGaussPointVecP() const { return &diff_gauss_points_vec; }

	vec3 GetPointVecByParam(double u0, double u1, double u2) const;
	vec3 GetPointVecByParam(double u0, double u1, double u2, double u3) const;

	bool ShareBasis(const Element&) const;
	bool IsVerticleElem() const;

	void SetMeshType(MeshType mesh_type) { this->mesh_type = mesh_type; }
	MeshType GetMeshType() const { return this->mesh_type; }
	void print(const string& str = "") const;
	string str() const;

};

inline int FindSubscriptOfFLayer(const vector<FLayer>& layer_list, const int idx)
{
	auto EqualFLayer = [idx](const FLayer& layer) {return layer.GetLayerIndex() == idx; }
	auto iter = find_if(layer_list.begin(), layer_list.end(), EqualFLayer);
	size_t diff = iter - layer_list.begin();
	return diff != layer_list.size() ? diff : -1;
}