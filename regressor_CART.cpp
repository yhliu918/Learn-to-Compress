// test accuracy of CART model
#include "../headers/common.h"
#include "../headers/file_manage.h"
#include "../headers/regress_tree.h"
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

	std::ofstream outfile("../reg_model_rle.txt", std::ios::out);
	read_csv("../train_data_try/training_rle_sample.txt", features, labels, 7933, 8);
	cout << "Start training" << endl;
	RegressionTree model;

	
	model.fit(features, labels);
	model.printTree();
	model.save(outfile);
	double accuracy = evaluate_model(model, features, labels, 10);
	cout << "Model MSE : " << accuracy << endl;
	
	/*
	model.rebuild(infile,4);
	double test_acc = model.test_pred( features_test, labels_test);
	cout << "test accuracy : " << test_acc << endl;
	*/
	
	return 0;
}
