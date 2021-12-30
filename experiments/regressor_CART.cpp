// test accuracy of CART model
#include "../headers/common.h"
#include "../headers/file_manage.h"
#include "../headers/regress_tree.h"
#include "../headers/model_selection.h"
#include "../headers/caltime.h"
using namespace std;
using namespace Eigen;
using namespace Codecset;

int main()
{
	MatrixXd features;
	VectorXd labels;
	MatrixXd features_test;
	VectorXd labels_test;

	//std::ifstream infile("../model.txt", std::ios::in);
	
	std::ofstream outfile("../reg_model/reg_model_piecewise.txt", std::ios::out);
	read_csv("../regression_model_train_data/training_piecewise.txt", features, labels, 5454, 8);
	cout << "Start training" << endl;
	RegressionTree model;

	double start = getNow();
	model.fit(features, labels);
	double end = getNow();
	model.printTree();
	model.save(outfile);
	double accuracy = evaluate_model(model, features, labels, 10);
	cout << "Model MSE : " << accuracy << endl;
	cout << "time: " << std::setprecision(8)
     << (end - start) << "s" << std::endl;
	
	/*
	model.rebuild(infile,4);
	double test_acc = model.test_pred( features_test, labels_test);
	cout << "test accuracy : " << test_acc << endl;
	*/
	
	return 0;
}
