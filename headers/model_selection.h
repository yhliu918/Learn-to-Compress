#pragma once
#include <vector>
#include <ctime>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;

template<class T>
double evaluate_model(T& model, const MatrixXd& df, const VectorXd& labels, int n_fold);

void split_to_folds(const MatrixXd& df, int n_fold, vector<vector<int>>& folds);

int unique_random(const vector<int>& unique, int range);

MatrixXd train_feature(const MatrixXd& df, const vector<vector<int>>& folds, int except);

VectorXd train_label(const VectorXd& labels, const vector<vector<int>>& folds, int except);

MatrixXd test_feature(const MatrixXd& df, const vector<vector<int>>& folds, int include);

VectorXd test_label(const VectorXd& labels, const vector<vector<int>>& folds, int include);

double calc_accuracy(const VectorXd& actual, const VectorXd& predicts);


template<class T>
double evaluate_model(T& model, const MatrixXd& df, const VectorXd& labels, int n_fold)
{
	vector<vector<int>> folds;
	split_to_folds(df, n_fold, folds);

	double accuracy = 0;
	for (int i = 0; i < n_fold; i++)
	{
		model.fit(train_feature(df, folds, i), train_label(labels, folds, i));

		VectorXd predicts = model.predict(test_feature(df, folds, i));
		accuracy += calc_accuracy(test_label(labels, folds, i), predicts);
	}
	return accuracy / n_fold;
}

void split_to_folds(const MatrixXd& df, int n_fold, vector<vector<int>>& folds)
// folds have n_fold vectors and each vector contains the indices of rows
{
	srand((unsigned)time(NULL));

	folds.resize(n_fold, vector<int>());

	vector<int> unique;
	size_t fold_size = df.rows() / n_fold;

	for (int i = 0; i < n_fold; i++)
		for (int j = 0; j < fold_size; j++)
			folds[i].push_back(unique_random(unique, (int)df.rows()));
}

int unique_random(const vector<int>& unique, int range)
{
	bool isOverlap;
	int num;
	do
	{
		num = rand() % range;
		isOverlap = false;
		for (int i = 0; i < unique.size(); i++)
			if (unique[i] == num)
			{
				isOverlap = true;
				break;
			}
	} while (isOverlap);
	return num;
}

MatrixXd train_feature(const MatrixXd& df, const vector<vector<int>>& folds, int except)
{
	size_t size = (folds.size() - 1) * folds[0].size();
	MatrixXd feature(size, df.cols());

	int i = 0;
	for (int j = 0; j < folds.size(); j++)
	{
		if (j == except)
			continue;

		for (const int& idx : folds[j])
		{
			feature.row(i) = df.row(idx);
			i++;
		}
	}
	return feature;
}

VectorXd train_label(const VectorXd& labels, const vector<vector<int>>& folds, int except)
{
	size_t size = (folds.size() - 1) * folds[0].size();
	VectorXd train_labels(size, 1);

	int i = 0;
	for (int j = 0; j < folds.size(); j++)
	{
		if (j == except)
			continue;

		for (const int& idx : folds[j])
		{
			train_labels[i] = labels[idx];
			i++;
		}
	}
	return train_labels;
}

MatrixXd test_feature(const MatrixXd& df, const vector<vector<int>>& folds, int include)
{
	size_t size = folds[0].size();
	MatrixXd feature(size, df.cols());

	int i = 0;
	for (const int& idx : folds[include])
	{
		feature.row(i) = df.row(idx);
		i++;
	}
	return feature;
}

VectorXd test_label(const VectorXd& labels, const vector<vector<int>>& folds, int include)
{
	size_t size = folds[0].size();
	VectorXd test_labels(size, 1);

	int i = 0;
	for (const int& idx : folds[include])
	{
		test_labels[i] = labels[idx];
		i++;
	}
	return test_labels;
}

double calc_accuracy(const VectorXd& actual, const VectorXd& predicts)
{
	double correct = 0;
	for (int i = 0; i < actual.size(); i++)
		if (actual[i] == predicts[i])
			correct++;
	return correct / actual.size();
}