
#ifndef PIECEWISEFIX_STRING_H_
#define PIECEWISEFIX_STRING_H_

#include "../common.h"
#include "../codecs.h"
#include "../bit_read.h"
#include "../bit_write.h"
#include "lr_string.h"
#include "string_utils.h"
#include "bit_read_string.h"
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/cpp_int/serialize.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/stream.hpp>
#define INF 0x7f7fffff

namespace Codecset
{
    // long_int to_bigint(uint8_t *v)
    // {
    //     namespace io = boost::iostreams;
    //     namespace ba = boost::archive;

    //     long_int i;
    //     {
    //         std::vector<char> chars{v, v+sizeof(i)};
    //         io::stream_buffer<io::array_source> bb(chars.data(), chars.size());
    //         ba::binary_iarchive ia(bb, ba::no_header | ba::no_tracking | ba::no_codecvt);
    //         ia >> i;
    //     }
    //     return i;
    // }
    // bit + theta0 + theta1 + delta

    uint8_t *encodeArray8_string(std::vector<std::string> &string_vec, const size_t length, uint8_t *res, size_t nvalue/*, long_int & theta0, long_int & theta1*/)
    {
        uint8_t *out = res;
        std::vector<long_int> ascii_vec;
        std::vector<long_int> index;
        for (int i = 0; i < string_vec.size(); i++)
        {
            ascii_vec.emplace_back(convertToASCII(string_vec[i]));
            index.emplace_back(long_int(i));
        }

        string_lr mylr;
        mylr.caltheta(index, ascii_vec, length);

        long_float theta0 = mylr.theta0;
        long_int theta0_int = theta0.convert_to<long_int>();
        long_float theta1 = mylr.theta1;
        long_int theta1_int = theta1.convert_to<long_int>();
        //std::cout << theta0_int << " " << theta1_int << std::endl;

        std::vector<long_int> delta;
        long_int max_delta = 0;
        for (auto i = 0; i < string_vec.size(); i++)
        {
            delta.emplace_back(ascii_vec[i] - (theta1_int * i + theta0_int));
            if (abs(delta[i]) > max_delta)
            {
                max_delta = abs(delta[i]);
            }
        }
        
        uint32_t max_bit = 0;
        if (max_delta)
        {
            max_bit = bits_long(max_delta) + 1;
        }
        //std::cout << "max bit " << max_bit << std::endl;
        memcpy(out, &max_bit, sizeof(uint32_t));
        out += sizeof(uint32_t);
        //models[nvalue*2] = theta0;
        //models[nvalue*2+1] = theta1;
        mpz_t z;
        mpz_init(z);
        mpz_set(z, theta0_int.backend().data());
        mpz_export(out, NULL, -1, 1, 0, 0, z);
        out += sizeof(long_int);


        mpz_set(z, theta1_int.backend().data());
        mpz_export(out, NULL, -1, 1, 0, 0, z);
        out += sizeof(long_int);
        mpz_clear(z);

        out = write_string_delta_string(delta.data(), out, max_bit, length);

        return out;
    }

    void decodeArray8_string(uint8_t *in, const size_t length, long_int *out, size_t nvalue, std::vector<std::string>& string_vec)
    {

        uint32_t maxbits;
        uint8_t *tmpin = in;
        memcpy(&maxbits, tmpin, sizeof(uint32_t));
        tmpin+=sizeof(uint32_t);
        //theta0 = reinterpret_cast<long_int*>(tmpin);
        mpz_t tmp;
        mpz_init(tmp);
        mpz_import(tmp, sizeof(long_int), -1, sizeof(tmpin[0]), 0, 0, tmpin);
        long_int theta0(tmp);
        std::cout<<theta0<<std::endl;

        tmpin += sizeof(long_int);

        mpz_import(tmp, sizeof(long_int), -1, sizeof(tmpin[0]), 0, 0, tmpin);
        long_int theta1(tmp);
        std::cout<<theta1<<std::endl;
        mpz_clear(tmp);
        tmpin += sizeof(long_int);

        read_all_bit_fix_string(tmpin ,0,0, length, maxbits,theta1,theta0, out);
        for (int i = 0; i < length; i++)
        {
            std::string tmp_string = convertToString(out[0]);
            string_vec.emplace_back(tmp_string);
            out++;
        }

    }

} // namespace FastPFor

#endif /* SIMDFASTPFOR_H_ */
