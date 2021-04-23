#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;

// declaration

void read_csv(string filename, MatrixXd& df, VectorXd& labels, int nRow, int nCol);

void split_line(string& line, MatrixXd& df, int row, string& label);

void make_class_vector(const vector<string>& labelStr, VectorXd& labels);

int string_label_to_int(const string& label, const vector<string>& labelUniq);

// implementation

void read_csv(string filename, MatrixXd& df, VectorXd& labels, int nRow, int nCol)
// Warning: The last column is always a class vector
{
	ifstream fin(filename);

	if (!fin.is_open())
		cout << "Error(read_csv(string, MatrixXd&, VectorXd&, int, int)): File not found." << endl;
	else
	{
		df.resize(nRow, nCol - 1);

		string line, label;
		vector<string> labelStr;

		int i = 0;
		while (getline(fin, line))
		{
			if (line == "")
				continue;

			split_line(line, df, i, label);

			if (label == "/0" || i >= df.rows())
			{
				cout << "Error(read_csv(string, MatrixXd&, VectorXd&, int, int)): " <<
					"Something went wrong." << endl;
				break;
			}

			labelStr.push_back(label);
			i++;
		}

		labels.resize(nRow);
		make_class_vector(labelStr, labels);

		fin.close();
	}
}

void split_line(string& line, MatrixXd& df, int row, string& label)
{
	string word = "";
	int col = 0;
	for (const char& ch : line)
	{
		if (ch == ',')
		{
			df(row, col) = std::stod(word);
			word = "";
			col++;
		}
		else
			word += ch;
	}

	if (col != df.cols())
	{
		cout << "Error(split_line()): Missing value or too many values." << endl;
		label = "\0";
	}
	else
		label = word;
}

void make_class_vector(const vector<string>& labelStr, VectorXd& labels)
{
	vector<string> labelUniq;
	int label;
	for (int i = 0; i < labels.size(); i++)
	{
		label = string_label_to_int(labelStr[i], labelUniq);

		if (label == -1)
		{
			labels[i] = (int)labelUniq.size();
			labelUniq.push_back(labelStr[i]);
		}
		else
			labels[i] = label;
	}
}

int string_label_to_int(const string& label, const vector<string>& labelUniq)
{
	unsigned i;
	for (i = 0; i < labelUniq.size(); i++)
		if (labelUniq[i] == label)
			return i;
	return -1;
}