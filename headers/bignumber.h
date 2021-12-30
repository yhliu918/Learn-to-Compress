#ifndef BigNumber_H
#define BigNumber_H

#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
constexpr auto precision = 50;


class BigNumber
{
public:
	template<class T>BigNumber(const T& value);
	
	bool isNegative() const;
	bool isPositive() const;
	bool isDecimal() const;
	BigNumber negate();
	BigNumber trimLeadingZeros();
	BigNumber normalize();
	int getPointIndex();
	// length of integer
	int getIntLength();
	// length of decimal (including decimal point)
	int getDecLength();
	std::string getString();
	BigNumber setString(const std::string& newStr);
	

	BigNumber add(BigNumber other);
	BigNumber subtract(BigNumber other);
	BigNumber multiply(BigNumber other);
	BigNumber divide(BigNumber other);
	bool equals(const BigNumber& other);

	friend BigNumber operator+(BigNumber bn1, const BigNumber& bn2);
	friend BigNumber operator-(BigNumber bn1, const BigNumber& bn2);
	friend BigNumber operator*(BigNumber bn1, const BigNumber& bn2);
	friend BigNumber operator/(BigNumber bn1, const BigNumber& bn2);
	friend bool operator==(BigNumber bn1, const BigNumber& bn2);

	friend std::ostream& operator<<(std::ostream& os, const BigNumber& num);
	friend bool operator>(BigNumber bn1, const BigNumber& bn2);
	friend bool operator<(BigNumber bn1, const BigNumber& bn2);
	friend bool operator>=(BigNumber bn1, const BigNumber& bn2);
	friend bool operator<=(BigNumber bn1, const BigNumber& bn2);

	BigNumber& operator=(const BigNumber& other);
	BigNumber& operator+=(const BigNumber& other);
	BigNumber& operator-=(const BigNumber& other);
	BigNumber& operator*=(const BigNumber& other);
	BigNumber& operator/=(const BigNumber& other);

	BigNumber& operator++();
	BigNumber& operator--();

	unsigned int operator[](int index);
private:
	std::string _numberString;
};


template<class T>
inline BigNumber::BigNumber(const T& value)
{
	std::stringstream ss;
	ss << value;
	std::string sout;
	sout = ss.str();
	if (sout == "") {
		sout = "0";
		/*std::cerr << "Warning: Initialization variable cannot be blank !" << std::endl;
		return;*/
	}
	if (sout.find_first_not_of('0') != std::string::npos) {
		sout.assign(sout.begin() + sout.find_first_not_of('0'), sout.end());
	}
	if (sout.find_first_not_of('.') == 1) {
		sout.insert(sout.begin(), '0');
	}
	if (sout.find_last_not_of('0') != std::string::npos) {
		if (sout.find_first_of('.') != std::string::npos) {
			if (sout.find_first_not_of('.') != std::string::npos)
				sout.assign(sout.begin(), sout.begin() + sout.find_last_not_of('0') + 1);
		}
	}
	if (sout.at(sout.size() - 1) == '.') {
		sout.erase(sout.size() - 1, sout.size());
	}
	if (sout.find("-.") != std::string::npos) {
		sout.insert(sout.find("-.") + 1, "0");
	}
	this->_numberString = sout;
}

#endif