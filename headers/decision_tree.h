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

struct Node
{
	Question Q;
	Node* left = nullptr;
	Node* right = nullptr;
	vector<int> labels;
	~Node()
	{
		delete left;
		delete right;
	}
};

class DecisionTree
{
private:
	Node* root;
public:
	DecisionTree();
	~DecisionTree();
	void fit(const MatrixXd& X, const VectorXd& Y);
	VectorXd predict(const MatrixXd& X);
	double test_pred(const MatrixXd& X, const VectorXd& Y);
	void printTree();
	void save(ofstream &ff);
	void rebuild(ifstream &ff,int classes);
};

namespace dt
{
	vector<int> init_split(const VectorXd& Y);

	vector<int> init_cols(int size);

	void build_tree(Node*& node, const MatrixXd& X, const VectorXd& Y, const vector<int>& split,
		const vector<int>& cols);

	Question find_best_question(const MatrixXd& X, const VectorXd& Y, const vector<int>& split,
		const vector<int>& cols, double& best_gain);

	double gini(const VectorXd& Y, const vector<int>& split);

	vector<int> count_class(const VectorXd& Y, const vector<int>& split);

	int n_classes(const VectorXd& Y);

	vector<double> unique_values(const VectorXd& col, const vector<int>& split);

	bool isOverlap(const vector<double>& unique, double value);

	vector<vector<int>> split_node(const Question& Q, const MatrixXd& X, const vector<int>& split);

	double info_gain(const VectorXd& Y, const vector<int>& left, const vector<int>& right, double current);

	vector<int> erase_taken_col(const Question& Q, const vector<int>& cols);

	double predict_implementation(const RowVectorXd& x, Node* node);

	void print_implementation(Node* node, int64_t width);

	void rebuild_implement(std::ifstream & ff, Node*& node,int class_num);

	void save_implement(std::ofstream & ff ,Node*& node);
}

DecisionTree::DecisionTree() : root(nullptr) {}

DecisionTree::~DecisionTree() { delete root; }

void DecisionTree::fit(const MatrixXd& X, const VectorXd& Y)
{
	dt::build_tree(root, X, Y, dt::init_split(Y), dt::init_cols((int)X.cols()));
}

vector<int> dt::init_split(const VectorXd& Y)
{
	vector<int> split;
	for (int i = 0; i < Y.size(); i++)
		split.push_back(i);
	return split;
}

vector<int> dt::init_cols(int size)
{
	vector<int> cols;
	for (int i = 0; i < size; i++)
		cols.push_back(i);
	return cols;
}

void dt::build_tree(Node*& node, const MatrixXd& X, const VectorXd& Y, const vector<int>& split,
	const vector<int>& cols)
{
	double gain;
	Question Q = find_best_question(X, Y, split, cols, gain);

	node = new Node;
	node->labels = count_class(Y, split);

	if (gain >= 0.08)
	{
		node->Q = Q;

		vector<vector<int>> splits = split_node(Q, X, split);
		vector<int> new_cols = erase_taken_col(Q, cols);

		build_tree(node->left, X, Y, splits[0], new_cols);
		build_tree(node->right, X, Y, splits[1], new_cols);
	}
}

Question dt::find_best_question(const MatrixXd& X, const VectorXd& Y, const vector<int>& split,
	const vector<int>& cols, double& best_gain)
{
	Question best_Q;
	best_gain = 0;

	double current_uncertainty = gini(Y, split);
	for (int idx : cols)
	{
		vector<double> unique = unique_values(X.col(idx), split);
		for (double value : unique)
		{
			Question Q(idx, value);
			vector<vector<int>> splits = split_node(Q, X, split);

			if (splits[0].size() == 0 || splits[1].size() == 0)
				continue;

			double gain = info_gain(Y, splits[0], splits[1], current_uncertainty);

			if (gain >= best_gain)
			{
				best_gain = gain;
				best_Q = Q;
			}
		}
	}
	return best_Q;
}

double dt::gini(const VectorXd& Y, const vector<int>& split)
{
	vector<int> counts = count_class(Y, split);

	double impurity = 1;
	for (auto count : counts)
		impurity -= std::pow(count / (double)split.size(), 2);
	return impurity;
}

vector<int> dt::count_class(const VectorXd& Y, const vector<int>& split)
{
	vector<int> counts(n_classes(Y));
	for (int idx : split)
	{
		int label = (int)Y[idx];
		counts[label]++;
	}
	return counts;
}

int dt::n_classes(const VectorXd& Y)
{
	return int(*std::max_element(Y.data(), Y.data() + Y.size()) + 1);
}

vector<double> dt::unique_values(const VectorXd& col, const vector<int>& split)
{
	vector<double> unique;
	for (auto idx : split)
		if (!isOverlap(unique, col[idx]))
			unique.push_back(col[idx]);
	return unique;
}

bool dt::isOverlap(const vector<double>& unique, double value)
{
	auto iter = std::find(unique.begin(), unique.end(), value);
	return iter != unique.end();
}

vector<vector<int>> dt::split_node(const Question& Q, const MatrixXd& X, const vector<int>& split)
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

double dt::info_gain(const VectorXd& Y, const vector<int>& left, const vector<int>& right, double current)
{
	double P = (double)left.size() / (left.size() + right.size());
	return current - P * gini(Y, left) - (1 - P) * gini(Y, right);
}

vector<int> dt::erase_taken_col(const Question& Q, const vector<int>& cols)
{
	vector<int> new_cols;
	for (int idx : cols)
		if (idx != Q.getCol())
			new_cols.push_back(idx);
	return new_cols;
}

VectorXd DecisionTree::predict(const MatrixXd& X)
{
	VectorXd labels(X.rows());
	for (int i = 0; i < X.rows(); i++)
		labels[i] = dt::predict_implementation(X.row(i), root);
	return labels;
}

double DecisionTree::test_pred(const MatrixXd& X, const VectorXd& Y){

	VectorXd pred(X.rows());
	pred = predict(X);
	int total=0;
	for (int i = 0; i < X.rows(); i++){
		if(pred[i]==Y[i]){
				total++;

		}
	}
	return (double)total/X.rows();

}
double dt::predict_implementation(const RowVectorXd& x, Node* node)
{
	while (node->left != nullptr && node->right != nullptr)
	{
		if (node->Q.match(x))
			node = node->left;
		else
			node = node->right;
	}

	auto max = std::max_element(node->labels.begin(), node->labels.end());
	return (double)std::distance(node->labels.begin(), max);
}

void DecisionTree::printTree() { dt::print_implementation(root, 0); }

void DecisionTree::rebuild(std::ifstream & ff,int class_num){
	
	dt::rebuild_implement(ff,root,class_num);
	printTree();

}
void dt::rebuild_implement(std::ifstream & ff, Node*& node,int class_num){
	if(ff.eof()){
		return ;
	}
	char c;
	ff>>c;
	node = new Node;
	
	if (c=='Q'){
		
		int tmpcol;
		double tmpval;
		
		ff>>tmpcol>>tmpval;
		
		Question tmpq(tmpcol,tmpval);
		node->Q=tmpq;
		dt::rebuild_implement(ff,node->left,class_num);
		dt::rebuild_implement(ff,node->right,class_num);
	}
	if(c=='L'){
		
		int val;
		ff>>val;
		std::vector<int> labeltmp(class_num);
		for(int i=0;i<class_num;i++){
			labeltmp[i]=0;

		}
		labeltmp[val]=1;
		node->labels = labeltmp;
		return;

	}
	

}
void DecisionTree::save(std::ofstream & ff ){
	dt::save_implement(ff,root);

}
void dt::save_implement(std::ofstream & ff ,Node*& node){
	// Q attribution value / # / L label
	if (node->left == nullptr && node->right == nullptr){
		auto max = std::max_element(node->labels.begin(), node->labels.end());
		ff<<'L'<<"  "<<std::distance(node->labels.begin(), max)<<std::endl;
		return;

	}
	else{
		ff<<'Q'<<"  "<<node->Q.getCol()<<"  "<<node->Q.getValue()<<std::endl;
		dt::save_implement(ff,node->left);
		dt::save_implement(ff,node->right);

	}


}
void dt::print_implementation(Node* node, int64_t width)
{
	if (node->left == nullptr && node->right == nullptr)
	{
		cout << setw(width + 4) << " " << "Predict : {";
		for (int i = 0; i < (int)node->labels.size(); i++)
			cout << "'" << i << "' : " << node->labels[i] << ", ";
		cout << "}" << endl;
	}
	else
	{
		cout << setw(width) << " " << "Q : X" << node->Q.getCol() + 1 << " >= " <<
			node->Q.getValue() << " ? " << endl;

		cout << setw(width) << " " << "--> True: " << endl;
		print_implementation(node->left, width + 4);

		cout << setw(width) << " " << "--> False: " << endl;
		print_implementation(node->right, width + 4);
	}
}
