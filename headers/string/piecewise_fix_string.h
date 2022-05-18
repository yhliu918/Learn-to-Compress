
#ifndef PIECEWISEFIX_STRING_H_
#define PIECEWISEFIX_STRING_H_

#include "../common.h"
#include "../codecs.h"
#include "../bit_read.h"
#include "../bit_write.h"
#include "lr_string.h"
#include "string_utils.h"
#include "bit_read_string.h"
#define INF 0x7f7fffff

namespace Codecset
{

    uint8_t *encodeArray8_string(std::vector<std::string> &string_vec,int start_idx, const size_t length, uint8_t *res, size_t nvalue/*, long_int & theta0, long_int & theta1*/)
    {
        uint8_t *out = res;
        std::vector<long_int> ascii_vec;
        std::vector<int> index;
        for (int i = 0; i < length; i++)
        {
            ascii_vec.emplace_back(convertToLongInt(string_vec[i+start_idx]));
            index.emplace_back(i);
        }

        string_lr mylr;
        mylr.caltheta(index, ascii_vec, length);

        std::vector<long_int> delta;
        long_int max_delta = 0;
        for (auto i = 0; i < length; i++)
        {
            long_int tmp_val = ascii_vec[i] - (mylr.theta1 * i + mylr.theta0);

            delta.emplace_back(tmp_val);   
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
        memcpy(out, &max_bit, sizeof(uint32_t));
        out += sizeof(uint32_t);

        mpz_t z;
        mpz_init(z);
        mpz_set(z, mylr.theta0.backend().data());
        auto theta0_len = (mpz_sizeinbase(z, 2) + 7) / 8;
        memcpy(out, &theta0_len, sizeof(uint32_t));
        out += sizeof(uint32_t);
        mpz_export(out, &theta0_len, -1, 1, 0, 0, z);
        out += theta0_len;


        mpz_set(z, mylr.theta1.backend().data());
        auto theta1_len = (mpz_sizeinbase(z, 2) + 7) / 8;
        memcpy(out, &theta1_len, sizeof(uint32_t));
        out += sizeof(uint32_t);
        mpz_export(out, &theta1_len, -1, 1, 0, 0, z);
        out += theta1_len;
        mpz_clear(z);
        if (max_bit)
        {
            out = write_string_delta_string(delta.data(), out, max_bit, length);
        }

        return out;
    }

    
    uint8_t *encodeArray8_string_128(std::vector<std::string> &string_vec,int start_idx, const size_t length, uint8_t *res, size_t nvalue/*, long_int & theta0, long_int & theta1*/)
    {
        uint8_t *out = res;
        std::vector<uint128_t> ascii_vec;
        std::vector<long_int> long_int_vec;
        std::vector<int> index;
        for (int i = 0; i < length; i++)
        {
            ascii_vec.emplace_back(convertTo128(string_vec[i+start_idx]));
            long_int_vec.emplace_back(convertToLongInt(string_vec[i+start_idx]));
            index.emplace_back(i);
            
        }

        string_lr mylr;
        mylr.caltheta(index, long_int_vec, length);

        uint128_t theta0=0;
        // auto theta0_len = (mpz_sizeinbase(mylr.theta0.backend().data(), 2) + 7) / 8;
        mpz_export(&theta0, NULL, -1, 1, 1, 0, mylr.theta0.backend().data());
        // print_u128_u(theta0);
        // printf("\n");

        uint128_t theta1=0;
        // auto theta1_len = (mpz_sizeinbase(mylr.theta1.backend().data(), 2) + 7) / 8;
        mpz_export(&theta1, NULL, -1, 1, 1, 0, mylr.theta1.backend().data());
        // print_u128_u(theta1);
        // printf("\n");
        // std::cout<<"theta0: "<<mylr.theta0<<std::endl;
        // std::cout<<"theta1: "<<mylr.theta1<<std::endl;

        std::vector<int128_t> delta;
        uint128_t max_delta = 0;
        for (auto i = 0; i < length; i++)
        {
            int128_t tmp_val = ascii_vec[i] - (theta1 * i + theta0);
            
            delta.emplace_back(tmp_val);   
            if (tmp_val<0)
            {
                tmp_val = -tmp_val;
            }

            if (tmp_val> max_delta)
            {
                max_delta = tmp_val;
            }
        }
        
        uint32_t max_bit = 0;
        if (max_delta)
        {
            max_bit = bits_128(max_delta) + 1;
        }
        // std::cout<< "max_bit: " << max_bit << std::endl;
        memcpy(out, &max_bit, sizeof(uint32_t));
        out += sizeof(uint32_t);

        memcpy(out, &theta0, sizeof(uint128_t));
        out += sizeof(uint128_t);
        memcpy(out, &theta1, sizeof(uint128_t));
        out += sizeof(uint128_t);
        if (max_bit)
        {
            out = write_string_delta_string_128(delta.data(), out, max_bit, length);
        }

        return out;
    }


    
    void decodeArray8_string(uint8_t *in, const size_t length,  size_t nvalue, std::vector<std::string>& string_vec)
    {
        uint32_t maxbits;
        uint8_t *tmpin = in;
        memcpy(&maxbits, tmpin, sizeof(uint32_t));
        //std::cout<<"max bit "<<maxbits<<std::endl;
        tmpin+=sizeof(uint32_t);

        //double start = getNow();
        uint32_t theta0_len;
        memcpy(&theta0_len, tmpin, sizeof(uint32_t));
        tmpin+=sizeof(uint32_t);
        //theta0 = reinterpret_cast<long_int*>(tmpin);
        mpz_t tmp;
        mpz_init(tmp);
        mpz_import(tmp, theta0_len, -1, 1, 0, 0, tmpin);
        long_int theta0(tmp);
        //std::cout<<theta0<<std::endl;

        tmpin += theta0_len;


        uint32_t theta1_len;
        memcpy(&theta1_len, tmpin, sizeof(uint32_t));
        tmpin+=sizeof(uint32_t);

        mpz_import(tmp, theta1_len, -1, 1, 0, 0, tmpin);
        long_int theta1(tmp);
        
        mpz_clear(tmp);
        tmpin += theta1_len;
        if (maxbits == 0)
        {
            for (int i = 0; i < length; i++)
            {
                long_int tmp_val = theta0 + theta1 * i;
                //string_vec.emplace_back(convertToString(&tmp_val));
            }
        }
        else
        {
            read_all_bit_fix_string(tmpin ,0,0, length, maxbits,theta1,theta0, string_vec);
        }
        // std::cout << "read_all_bit_fix_string time per int: " << std::setprecision(8)
        //     << (end-start) / length * 1000000000 << "ns" << std::endl;

        // start = getNow();
        // for (int i = 0; i < length; i++)
        // {
        //     std::string tmp_string = convertToString(*out[i]);
        //     string_vec.emplace_back(tmp_string);
        // }
        // end = getNow();
        // std::cout << "ascii to string time per int: " << std::setprecision(8)
        //     << (end-start) / length * 1000000000 << "ns" << std::endl;

    }

    void randomdecodeArray8(uint8_t *in, int idx,std::string& result){
        uint8_t *tmpin = in;
        uint32_t maxbits;
        memcpy(&maxbits, tmpin, sizeof(uint32_t));
        tmpin+=sizeof(uint32_t);
        uint32_t theta0_len;
        memcpy(&theta0_len, tmpin, sizeof(uint32_t));
        tmpin+=sizeof(uint32_t);

        mpz_t tmp;
        mpz_init(tmp);
        mpz_import(tmp, theta0_len, -1, 1, 0, 0, tmpin);
        long_int theta0(tmp);
        tmpin += theta0_len;

        uint32_t theta1_len;
        memcpy(&theta1_len, tmpin, sizeof(uint32_t));
        tmpin+=sizeof(uint32_t);
        mpz_import(tmp, theta1_len, -1, 1, 0, 0, tmpin);
        long_int theta1(tmp);
        tmpin += theta1_len;

        
        result = read_bit_fix_string_(tmpin ,maxbits, idx, theta1,theta0);
        mpz_clear(tmp);
        
    }


    void randomdecodeArray8_longint(uint8_t *in, int idx, mpz_t * result){
        uint8_t *tmpin = in;
        uint32_t maxbits;
        memcpy(&maxbits, tmpin, sizeof(uint32_t));
        tmpin+=sizeof(uint32_t);
        uint32_t theta0_len;
        memcpy(&theta0_len, tmpin, sizeof(uint32_t));
        tmpin+=sizeof(uint32_t);

        mpz_t tmp;
        mpz_init(tmp);
        mpz_import(tmp, theta0_len, -1, 1, 0, 0, tmpin);
        long_int theta0(tmp);
        tmpin += theta0_len;

        uint32_t theta1_len;
        memcpy(&theta1_len, tmpin, sizeof(uint32_t));
        tmpin+=sizeof(uint32_t);
        mpz_import(tmp, theta1_len, -1, 1, 0, 0, tmpin);
        long_int theta1(tmp);
        tmpin += theta1_len;

        read_bit_fix_string_long_int(tmpin ,maxbits, idx, theta1,theta0, result);
        //mpz_clear(tmp);

        // mpz_clear(tmp);
        
    }

     void randomdecodeArray8_longint_128(uint8_t *in, int idx, uint128_t * result){
        uint8_t *tmpin = in;
        uint32_t maxbits;
        memcpy(&maxbits, tmpin, sizeof(uint32_t));
        tmpin+=sizeof(uint32_t);

        uint128_t theta0;
        memcpy(&theta0, tmpin, sizeof(uint128_t));
        tmpin+=sizeof(uint128_t);

        uint128_t theta1;
        memcpy(&theta1, tmpin, sizeof(uint128_t));
        tmpin+=sizeof(uint128_t);
        read_bit_fix_string_long_int_128(tmpin ,maxbits, idx, theta1,theta0, result);
        //mpz_clear(tmp);

        // mpz_clear(tmp);
        
    }

} // namespace FastPFor

#endif /* SIMDFASTPFOR_H_ */
