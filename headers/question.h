#pragma once
#include <Eigen/Dense>
using namespace Eigen;

class Question
{
private:
	int col;
	double value;
public:
	Question();
	Question(int col, double value);
	bool match(const RowVectorXd& x) const;
	int getCol() const;
	double getValue() const;
	Question& operator=(const Question& other);
};

Question::Question() : col(0), value(0) {}

Question::Question(int col, double value) : col(col), value(value) {}

bool Question::match(const RowVectorXd& x) const
{
	return x[col] >= value;
}

int Question::getCol() const
{
	return col;
}

double Question::getValue() const
{
	return value;
}

Question& Question::operator=(const Question& other)
{
	col = other.col;
	value = other.value;
	return *this;
}