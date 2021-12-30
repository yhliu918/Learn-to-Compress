#include "../headers/bignumber.h"

#include <math.h>


std::ostream& operator<<(std::ostream& os, const BigNumber& num) {
	os << num._numberString;
	return os;
}

BigNumber operator+(BigNumber b1, const BigNumber& b2) {
    return b1.add(b2);
}

BigNumber operator-(BigNumber b1, const BigNumber& b2)
{
    return b1.subtract(b2);
}

BigNumber operator*(BigNumber bn1, const BigNumber& bn2)
{
    return bn1.multiply(bn2);
}

BigNumber operator/(BigNumber bn1, const BigNumber& bn2)
{
    return bn1.divide(bn2);
}

bool operator==(BigNumber b1, const BigNumber& b2) {
    return b1.equals(b2);
}

bool operator>(BigNumber bn1, const BigNumber& bn2)
{
    if (bn1.isNegative() || bn2.isNegative()) {
        if (bn1.isNegative() && bn2.isNegative()) {
            BigNumber bt = bn2;
            bn1._numberString.erase(0, 1);
            bt._numberString.erase(0, 1);
            return bn1 < bt;
        }
        else {
            return !(bn1.isNegative() && !bn2.isNegative());
        }
    }
    BigNumber bn3(bn2);
    int length_int1 = bn1.getIntLength();
    int length_int2 = bn3.getIntLength();

    int length_dec1 = bn1.getDecLength();
    int length_dec2 = bn3.getDecLength();

    // compare integer bits
    if (length_int1 >= length_int2) {
        if (length_int1 == length_int2) {
            for (int i = 0; i < length_int1; i++) {
                if (bn1._numberString[i] > bn3._numberString[i])
                    return true;
                else if (bn1._numberString[i] == bn3._numberString[i]) 
                    continue;
                else
                    return false;
            }
        }
        else
            return true;
    }
    else {
        return false;
    }

    // if integer parts equal, compare decimal parts
    int max_dec_length = (length_dec1 > length_dec2) ? length_dec1: length_dec2;
    for (int i = 1; i < max_dec_length; ++i) {
        // avoid exceeding the string range
        if (length_int1 + i > bn1._numberString.size())
            return false;
        if (length_int2 + i > bn3._numberString.size())
            return true;

        if (bn1._numberString[length_int1 + i] > bn3._numberString[length_int2 + i]) {
            return true;
        }
        else if (bn1._numberString[length_int1 + i] == bn3._numberString[length_int2 + i]) {
            continue;
        }
        else {
            return false;
        }
    }
    // bn1 = bn2
    return false;
}

bool operator<(BigNumber b1, const BigNumber& b2) {
    return !(b1 == b2) && !(b1 > b2);
}

bool operator>=(BigNumber bn1, const BigNumber& bn2)
{
    return bn1 > bn2 || bn1 == bn2;
}

bool operator<=(BigNumber bn1, const BigNumber& bn2)
{
    return bn1 < bn2 || bn1 == bn2;
}

BigNumber& BigNumber::operator=(const BigNumber& other)
{
    this->_numberString = other._numberString;
    return *this;
}

BigNumber& BigNumber::operator+=(const BigNumber& other)
{
    *this = *this + other;
    return *this;
}

BigNumber& BigNumber::operator-=(const BigNumber& other)
{
    *this = (*this).normalize() - other;
    return *this;
}

BigNumber& BigNumber::operator*=(const BigNumber& other)
{
    *this = *this * other;
    return *this;
}

BigNumber& BigNumber::operator/=(const BigNumber& other)
{
    *this = *this / other;
    return *this;
}

BigNumber& BigNumber::operator++()
{
    *this += BigNumber("1");
    return *this;
}

BigNumber& BigNumber::operator--()
{
    *this -= BigNumber("1");
    return *this;
}

unsigned int BigNumber::operator[](int index)
{
    if (this->_numberString[index] == '-') {
        std::cerr << "You cannot get the negative sign from the number" << std::endl;
    }
    return static_cast<unsigned int>(this->_numberString[index] - '0');
}

BigNumber BigNumber::add(BigNumber other)
{

    // make greater = bn1, less = bn2;
    BigNumber bn1 = other > * this ? other : *this;
    BigNumber bn2 = other > * this ? *this : other;

    if (bn1.isNegative() || bn2.isNegative()) {
        if (bn1.isNegative() && bn2.isNegative()) {
            return bn1.negate().add(bn2.negate()).negate();
        }
        else if (bn1.isNegative() && !bn2.isNegative()) {
            return bn1.negate().subtract(bn2).negate();
        }
        else {
            return bn2.negate().subtract(bn1).negate();
        }
    }

    std::string results;

    int diff_int = bn1.getIntLength() - bn2.getIntLength();
    // padding integer bits
    for (int i = 0; i < diff_int; ++i) {
        bn2._numberString.insert(bn2._numberString.begin(), '0');
    }
    // padding decimal bits
    if (bn1.getDecLength() > bn2.getDecLength()) {
        int diff_dec = bn1.getDecLength() - bn2.getDecLength();
        for (int i = 0; i < diff_dec; ++i) {
            if (i == 0 && !bn2.isDecimal())
                bn2._numberString.insert(bn2._numberString.end(), '.');
            else
                bn2._numberString.insert(bn2._numberString.end(), '0');
        }
    }
    else {
        int diff_dec = bn2.getDecLength() - bn1.getDecLength();
        for (int i = 0; i < diff_dec; ++i) {
            if (i == 0 && !bn1.isDecimal())
                bn1._numberString.insert(bn1._numberString.end(), '.');
            else
                bn1._numberString.insert(bn1._numberString.end(), '0');
        }
    }
    if (bn1.isDecimal() || bn2.isDecimal()) results.insert(results.end(), '.');

    int min_dec_length = bn1.getIntLength(); // the same as bn2 integer length
    int max_dec_length = bn1._numberString.size(); // the same as bn2 string length
    int carry = 0;
    // calculate decimal bits
    for (int i = max_dec_length - 1; i > min_dec_length; --i) {
        int sum = (bn1._numberString[i] - '0') + (bn2._numberString[i] - '0') + carry;
        carry = 0;
        if (sum <= 9 || i == 0) {
            results.insert(1, std::to_string(sum));
        }
        else {
            results.insert(1, std::to_string(sum % 10));
            carry = 1;
        }
    }
    int max_int_length = bn1.getIntLength();
    // calculate integer bits
    for (int i = max_int_length - 1; i >= 0; --i) {
        int sum = (bn1._numberString[i] - '0') + (bn2._numberString[i] - '0') + carry;
        carry = 0;
        if (sum <= 9 || i == 0) {
            results.insert(0, std::to_string(sum));
        }
        else {
            results.insert(0, std::to_string(sum % 10));
            carry = 1;
        }
    }
    
    return BigNumber(results).normalize();
}

BigNumber BigNumber::subtract(BigNumber other)
{
    BigNumber bn1 = *this, bn2 = other;
    if (bn1.isNegative() || bn2.isNegative()) {
        if (bn1.isNegative() && bn2.isNegative()) {
            return bn1.negate().add(bn2.negate()).negate().normalize();
        }
        else if (bn1.isNegative() && !bn2.isNegative()) {
            return bn1.negate().add(bn2).negate().normalize();
        }
        else {
            return bn2.negate().add(bn1).normalize();
        }
    }
  
    if (bn1 < bn2) {
        //Negative answer
        std::string t = bn2.subtract(*this).negate().getString();
       /* for (unsigned int i = 1; i < t.length(); ++i) {
            if (t[i] != '0') break;
            t.erase(1, 1);
        }*/
        return BigNumber(t);
    }

    //This next if-block fixes the case where the digit difference is greater than 1
    //100 - 5 is an example. This code adds 0's to make it, for example, 100 - 05, which
    //allows the rest of the subtraction code to work.

    int diff_int = bn1.getIntLength() - bn2.getIntLength();
    // padding integer bits
    for (int i = 0; i < diff_int; ++i) {
        bn2._numberString.insert(bn2._numberString.begin(), '0');
    }
    // padding decimal bits
    if (bn1.getDecLength() > bn2.getDecLength()) {
        int diff_dec = bn1.getDecLength() - bn2.getDecLength();
        for (int i = 0; i < diff_dec; ++i) {
            if (i == 0 && !bn2.isDecimal())
                bn2._numberString.insert(bn2._numberString.end(), '.');
            else
                bn2._numberString.insert(bn2._numberString.end(), '0');
        }
    }
    else {
        int diff_dec = bn2.getDecLength() - bn1.getDecLength();
        for (int i = 0; i < diff_dec; ++i) {
            if (i == 0 && !bn1.isDecimal())
                bn1._numberString.insert(bn1._numberString.end(), '.');
            else
                bn1._numberString.insert(bn1._numberString.end(), '0');
        }
    }

    std::string complement, one_tmp;
    int n = 0, p = 0;
    for (int i = 0; i < (int)bn2._numberString.size(); i++) {
        if (bn2._numberString[i] == '.'){
            complement.insert(complement.size(), ".");
            one_tmp.insert(one_tmp.size(), ".");
        }
        else {
            switch (bn2._numberString[i])
            {
            case '0': n = 9; break;
            case '1': n = 8; break;
            case '2': n = 7; break;
            case '3': n = 6; break;
            case '4': n = 5; break;
            case '5': n = 4; break;
            case '6': n = 3; break;
            case '7': n = 2; break;
            case '8': n = 1; break;
            case '9': n = 0; break;
            }
            complement.insert(complement.size(), std::to_string(n));
            one_tmp.insert(one_tmp.size(), "0");
        }
    }
    one_tmp[one_tmp.size() - 1] = '1';
    BigNumber a(complement), b(one_tmp);
    BigNumber c = a + b;
    BigNumber d = bn1 + c;
    if (bn1.getIntLength() != d.getIntLength())
        d._numberString.erase(0, 1);
    return d.normalize();
}

BigNumber BigNumber::multiply(BigNumber other)
{
    BigNumber bn1 = other > * this ? other : *this;
    BigNumber bn2 = other > * this ? *this : other;
    if (bn1.isNegative() || bn2.isNegative()) {
        if (bn1.isNegative() && bn2.isNegative()) {
            return bn1.negate().multiply(bn2.negate());
        }
        else if (bn1.isNegative() && !bn2.isNegative()) {
            return bn1.negate().multiply(bn2).negate();
        }
        else {
            return bn2.negate().multiply(bn1).negate();
        }
    }
    if (bn1 == 0 || bn2 == 0) return 0;

    std::string calc_tmp;
    BigNumber results = "";
    int carry = 0;
    int mult_times = bn2._numberString.size() - 1; 
   
    int compat = 0;
    for (int i = mult_times; i >= 0; i--) {
        if (bn2._numberString[i] == '.' || bn2._numberString[i] == '0') {
            if (bn2._numberString[i] == '0')
                compat++;
            continue;
        }
            
        for (int d = (int)bn1._numberString.size() - 1; d >= 0; d--) {
            if (bn1._numberString[d] != '.') {
                int mult = (bn1._numberString[d] - '0') * (bn2._numberString[i] - '0') + carry;
                carry = 0;
                if (mult <= 9 || d == 0) {
                    calc_tmp.insert(0, std::to_string(mult));
                }
                else {
                    calc_tmp.insert(0, std::to_string(mult % 10));
                    carry = floor(mult/10);
                }
            }
        }
        for(int i=0;i< compat;i++)
            calc_tmp.append("0");
        results += calc_tmp;
        calc_tmp = "";
        compat++;
      
    }
    int a = (bn1.getDecLength() - 1) < 0 ? 0 : (bn1.getDecLength() - 1);
    int b = (bn2.getDecLength() - 1) < 0 ? 0 : (bn2.getDecLength() - 1);
    int c = results._numberString.size() - (a + b);
    if (c <  0)
        for (c; c <= 0; c++)
            results._numberString.insert(0, "0");
    results._numberString.insert(c, ".");
    
    return results.normalize();
}

BigNumber BigNumber::divide(BigNumber other)
{
    if (other == 0) {
        std::cerr << "You cannot divide by 0!" << std::endl;
    }
    BigNumber bn1 = *this, bn2 = other;
    bool sign = false;
    if (bn1.isNegative() && bn2.isNegative()) {
        bn1.negate();
        bn2.negate();
    }
    else if (bn1.isNegative() && !bn2.isNegative()) {
        bn1.negate();
        sign = true;
    }
    else if (!bn1.isNegative() && bn2.isNegative()) {
        bn2.negate();
        sign = true;
    }
    BigNumber quotient = 0;
    BigNumber dividend(bn1);
    while (bn1 >= bn2) {
        bn1 -= bn2;
        ++quotient;
    }
    if (sign) quotient.negate();
    if (quotient * bn2 == dividend)
        return quotient;
    BigNumber p = "0.1";
    BigNumber p_tmp = "0";
    BigNumber p_tmp2 = "0.1";
    BigNumber result = "0";
    BigNumber result_tmp = quotient;
    for (int i = 0; i < precision; i++) {
        for (int n = 0; n < 10; n++) {
            result = (result_tmp + p_tmp);
            if (result * bn2 <= dividend)
                p_tmp += p;
            else {
                result = result - p;
                break;
            }
        }
        p._numberString.insert(2, "0");
        p_tmp = p;
        result_tmp = result;
    }
   
    return result;
}

bool BigNumber::isNegative() const {
    return this->_numberString[0] == '-';
}
bool BigNumber::isPositive() const {
    return !this->isNegative();
}


bool BigNumber::isDecimal() const
{
    struct checkDecimal {
        bool operator()(const char& key) {
            return key == '.';
        };
    };
    std::string::const_iterator found = find_if(this->_numberString.begin(), this->_numberString.end(), checkDecimal());

    return (found != this->_numberString.end());
}

BigNumber BigNumber::negate() {
    if (this->_numberString[0] == '-') {
        this->_numberString.erase(0, 1);
    }
    else {
        this->_numberString.insert(this->_numberString.begin(), '-');
    }
    return *this;
}

BigNumber BigNumber::trimLeadingZeros() {
    BigNumber b = *this;
    if (b._numberString.find_first_not_of('0') != std::string::npos) {
        b.setString(b._numberString.erase(0, b._numberString.find_first_not_of('0')));
    }
    return b;
}

BigNumber BigNumber::normalize()
{
    BigNumber bn = *this;
    
    /*if (bn._numberString.find_first_not_of('0') != std::string::npos) {
        bn._numberString.assign(bn._numberString.begin() + bn._numberString.find_first_not_of('0'), bn._numberString.end());
    }
    if (bn._numberString.find_last_not_of('0') != std::string::npos) {
        if (bn._numberString.find_first_not_of('.') < bn._numberString.find_last_not_of('0'))
            bn._numberString.assign(bn._numberString.begin(), bn._numberString.begin() + bn._numberString.find_last_not_of('0') + 1);
    }
    if (bn._numberString.find_first_not_of('.') == 1) {
        bn._numberString.insert(bn._numberString.begin(), '0');
    }*/
    if (bn._numberString.find_first_not_of('0') != std::string::npos) {
        bn._numberString.assign(bn._numberString.begin() + bn._numberString.find_first_not_of('0'), bn._numberString.end());
    }
    if (bn._numberString.find_first_not_of('.') == 1) {
        bn._numberString.insert(bn._numberString.begin(), '0');
    }
    /*
    if (bn._numberString.find_last_not_of('0') != std::string::npos) {
        std::cout << bn._numberString.find_first_not_of('.') << std::endl;
        if (bn._numberString.find_first_not_of('.') != 0) {
            bn._numberString.assign(bn._numberString.begin(), bn._numberString.begin() + bn._numberString.find_last_not_of('0') + 1);
        }
        
        
    }*/

    if (bn._numberString.find_last_not_of('0') != std::string::npos) {
        if (bn._numberString.find_first_of('.') != std::string::npos) {
            if (bn._numberString.find_first_not_of('.') != std::string::npos)
                bn._numberString.assign(bn._numberString.begin(), bn._numberString.begin() + bn._numberString.find_last_not_of('0') + 1);
        }
        /*
        if (sout.find_first_not_of('.') < sout.find_last_not_of('0'))
            sout.assign(sout.begin(), sout.begin() + sout.find_last_not_of('0') + 1);*/
    }
    if (bn._numberString.at(bn._numberString.size() - 1) == '.') {
        bn._numberString.erase(bn._numberString.size() - 1, bn._numberString.size());
    }
    /*
    if (bn._numberString.find_last_not_of('.') != std::string::npos) {
        bn.setString(bn._numberString.erase(bn._numberString.find_last_not_of('.') + 1, bn._numberString.size()));
    }*/
    
    /*if (bn._numberString.find("-.")!= std::string::npos) {
        bn.setString(bn._numberString.insert(1, "0"));
    }*/
    return bn;
}

std::string BigNumber::getString()
{
    return this->_numberString;
}

BigNumber BigNumber::setString(const std::string& newStr)
{
    this->_numberString = newStr;
    return *this;
}

int BigNumber::getPointIndex()
{
    int index = 0;
    bool check = false;
    for (std::string::const_iterator key = this->_numberString.begin(); key != this->_numberString.end(); key++) {
        if ((*key) == '.') {
            check = true;
            break;
        }
        ++index;
    }
    return (check) ? index : -1;
}

int BigNumber::getIntLength()
{
    return (this->isDecimal()) ? this->getPointIndex() : this->_numberString.size();
}

int BigNumber::getDecLength()
{
    return this->_numberString.size() - this->getIntLength();
}

bool BigNumber::equals(const BigNumber& other)
{
    return this->_numberString == other._numberString;
}
