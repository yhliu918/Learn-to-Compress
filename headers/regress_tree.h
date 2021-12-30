#pragma once
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <Eigen/Dense>
#include<queue>
#include<fstream>
#include "question.h"
using namespace std;
using namespace Eigen;

struct Node_reg
{
	Question Q;
	Node_reg* left = nullptr;
	Node_reg* right = nullptr;
	double labels;
	~Node_reg()
	{
		delete left;
		delete right;
	}
};

class RegressionTree
{
private:
	Node_reg* root= nullptr;
public:
	RegressionTree();
	~RegressionTree();
	void fit(const MatrixXd& X, const VectorXd& Y);
	VectorXd predict(const MatrixXd& X);
	double test_pred(const MatrixXd& X, const VectorXd& Y);
	void printTree();
	void save(ofstream &ff);
	void rebuild(ifstream &ff,int classes);
};

namespace rt
{
	vector<int> init_split(const VectorXd& Y);

	vector<int> init_cols(int size);

	void build_tree(Node_reg*& node, const MatrixXd& X, const VectorXd& Y, const vector<int>& split,
		const vector<int>& cols);

	Question find_best_question(const MatrixXd& X, const VectorXd& Y, const vector<int>& split,
		const vector<int>& cols, double& best_gain);


	
	double mse(const VectorXd& Y, const vector<int>& split);

	double cal_avg(const VectorXd& Y, const vector<int>& split);

	int n_classes(const VectorXd& Y);

	vector<double> unique_values(const VectorXd& col, const vector<int>& split);

	bool isOverlap(const vector<double>& unique, double value);

	vector<vector<int>> split_Node_reg(const Question& Q, const MatrixXd& X, const vector<int>& split);

	double info_gain_reg(const VectorXd& Y, const vector<int>& left, const vector<int>& right, double current);

	vector<int> erase_taken_col(const Question& Q, const vector<int>& cols);

	double predict_implementation_reg(const RowVectorXd& x, Node_reg* node);

	void print_implementation(Node_reg* node, int64_t wirth);

	void rebuild_implement(std::ifstream & ff, Node_reg*& node,int class_num);

	void save_implement(std::ofstream & ff ,Node_reg*& node);
}

RegressionTree::RegressionTree() : root(nullptr) {}

RegressionTree::~RegressionTree() { //delete root; 
}

void RegressionTree::fit(const MatrixXd& X, const VectorXd& Y)
{
	rt::build_tree(root, X, Y, rt::init_split(Y), rt::init_cols((int)X.cols()));
}

vector<int> rt::init_split(const VectorXd& Y)
{
	vector<int> split;
	for (int i = 0; i < Y.size(); i++)
		split.push_back(i);
	return split;
}

vector<int> rt::init_cols(int size)
{
	vector<int> cols;
	for (int i = 0; i < size; i++)
		cols.push_back(i);
	return cols;
}

void rt::build_tree(Node_reg*& node, const MatrixXd& X, const VectorXd& Y, const vector<int>& split,
	const vector<int>& cols)
{
	double gain;
	Question Q = find_best_question(X, Y, split, cols, gain);

	node = new Node_reg;
	node->labels = cal_avg(Y, split);


	if (gain >= 0.08)
	{
		node->Q = Q;

		vector<vector<int>> splits = split_Node_reg(Q, X, split);
		vector<int> new_cols = erase_taken_col(Q, cols);

		build_tree(node->left, X, Y, splits[0], new_cols);
		build_tree(node->right, X, Y, splits[1], new_cols);
	}
}

Question rt::find_best_question(const MatrixXd& X, const VectorXd& Y, const vector<int>& split,
	const vector<int>& cols, double& best_gain)
{
	Question best_Q;
	best_gain = 0;

	double current_uncertainty = mse(Y, split);
	for (int idx : cols)
	{
		vector<double> unique = unique_values(X.col(idx), split);
		for (double value : unique)
		{
			Question Q(idx, value);
			vector<vector<int>> splits = split_Node_reg(Q, X, split);

			if (splits[0].size() == 0 || splits[1].size() == 0)
				continue;

			double gain = info_gain_reg(Y, splits[0], splits[1], current_uncertainty);

			if (gain >= best_gain)
			{
				best_gain = gain;
				best_Q = Q;
			}
		}
	}
	return best_Q;
}



double rt::mse(const VectorXd& Y, const vector<int>& split)
{

	double avg = cal_avg(Y,split);
	double sum = 0;
	for (int idx : split)
	{
		double truth = (double)Y[idx];
		sum += std::pow(truth-avg,2);
	}
	
	return sum;
}


double rt::cal_avg(const VectorXd& Y, const vector<int>& split)
{
	double sum = 0;
	for (int idx : split)
	{
		double truth = (double)Y[idx];
		sum+=truth;
	}
	double avg = sum/split.size();
	
	return avg;
}

int rt::n_classes(const VectorXd& Y)
{
	return int(*std::max_element(Y.data(), Y.data() + Y.size()) + 1);
}

vector<double> rt::unique_values(const VectorXd& col, const vector<int>& split)
{
	vector<double> unique;
	for (auto idx : split)
		if (!isOverlap(unique, col[idx]))
			unique.push_back(col[idx]);
	return unique;
}

bool rt::isOverlap(const vector<double>& unique, double value)
{
	auto iter = std::find(unique.begin(), unique.end(), value);
	return iter != unique.end();
}

vector<vector<int>> rt::split_Node_reg(const Question& Q, const MatrixXd& X, const vector<int>& split)
{
	vector<vector<int>> splits(2, vector<int>());
	for (int idx : split)
	{
		if (Q.match(X.row(idx)))
			splits[0].push_back(idx);
		else
			splits[1].push_back(idx);
	}
	return splits;
}


double rt::info_gain_reg(const VectorXd& Y, const vector<int>& left, const vector<int>& right, double current)
{
	return current - (mse(Y,left)+mse(Y,right));
}

vector<int> rt::erase_taken_col(const Question& Q, const vector<int>& cols)
{
	vector<int> new_cols;
	for (int idx : cols)
		if (idx != Q.getCol())
			new_cols.push_back(idx);
	return new_cols;
}

VectorXd RegressionTree::predict(const MatrixXd& X)
{
	VectorXd labels(X.rows());
	for (int i = 0; i < X.rows(); i++)
		labels[i] = rt::predict_implementation_reg(X.row(i), root);
	return labels;
}

double RegressionTree::test_pred(const MatrixXd& X, const VectorXd& Y){

	VectorXd pred(X.rows());
	pred = predict(X);
	double total=0;
	for (int i = 0; i < X.rows(); i++){
		total+=pow(pred[i]-Y[i],2);
	}
	return (double)total/X.rows();

}


double rt::predict_implementation_reg(const RowVectorXd& x, Node_reg* node)
{
	while (node->left != nullptr && node->right != nullptr)
	{
		if (node->Q.match(x))
			node = node->left;
		else
			node = node->right;
	}

	auto avg = node->labels;
	return avg;
}

void RegressionTree::printTree() { rt::print_implementation(root, 0); }

void RegressionTree::rebuild(std::ifstream & ff,int class_num){
	
	rt::rebuild_implement(ff,root,class_num);
	printTree();

}
void rt::rebuild_implement(std::ifstream & ff, Node_reg*& node,int class_num){
	if(ff.eof()){
		return ;
	}
	char c;
	ff>>c;
	node = new Node_reg;
	
	if (c=='Q'){
		
		int tmpcol;
		double tmpval;
		
		ff>>tmpcol>>tmpval;
		
		Question tmpq(tmpcol,tmpval);
		node->Q=tmpq;
		rt::rebuild_implement(ff,node->left,class_num);
		rt::rebuild_implement(ff,node->right,class_num);
	}
	if(c=='L'){
		
		double val;
		ff>>val;
		node->labels = val;
		return;

	}
	

}
void RegressionTree::save(std::ofstream & ff ){
	rt::save_implement(ff,root);

}
void rt::save_implement(std::ofstream & ff ,Node_reg*& node){
	// Q attribution value / # / L label
	if (node->left == nullptr && node->right == nullptr){
		auto max = node->labels;
		ff<<'L'<<"  "<<max<<std::endl;
		return;

	}
	else{
		ff<<'Q'<<"  "<<node->Q.getCol()<<"  "<<node->Q.getValue()<<std::endl;
		rt::save_implement(ff,node->left);
		rt::save_implement(ff,node->right);

	}


}
void rt::print_implementation(Node_reg* node, int64_t wirth)
{
	if (node->left == nullptr && node->right == nullptr)
	{
		cout << setw(wirth + 4) << " " << "Predict : {";
		cout << node->labels;
		cout << "}" << endl;
	}
	else
	{
		cout << setw(wirth) << " " << "Q : X" << node->Q.getCol() + 1 << " >= " <<
			node->Q.getValue() << " ? " << endl;

		cout << setw(wirth) << " " << "--> True: " << endl;
		print_implementation(node->left, wirth + 4);

		cout << setw(wirth) << " " << "--> False: " << endl;
		print_implementation(node->right, wirth + 4);
	}
}
