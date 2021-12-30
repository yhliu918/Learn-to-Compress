// test accuracy of CART model
#include "../headers/common.h"
#include "../headers/file_manage.h"
#include "../headers/decision_tree.h"
#include "../headers/model_selection.h"
using namespace std;
using namespace Eigen;

int main()
{
	MatrixXd features;
	VectorXd labels;
	MatrixXd features_test;
	VectorXd labels_test;

	//std::ifstream infile("../model.txt", std::ios::in);

	std::ofstream outfile("../model_6dim.txt", std::ios::out);
	read_csv("../train_data/training.csv", features, labels, 1663, 7);
	
	DecisionTree model;

	
	model.fit(features, labels);
	model.printTree();
	model.save(outfile);
	double accuracy = evaluate_model(model, features, labels, 30);
	cout << "Model accuracy : " << accuracy << endl;
	
	/*
	model.rebuild(infile,4);
	double test_acc = model.test_pred( features_test, labels_test);
	cout << "test accuracy : " << test_acc << endl;
	*/
	
	return 0;
}
