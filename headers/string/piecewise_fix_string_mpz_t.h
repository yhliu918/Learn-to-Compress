#ifndef PIECEWISEFIX_STRING_MPZ_H_
#define PIECEWISEFIX_STRING_MPZ_H_

#include "../common.h"
#include "../codecs.h"
#include "../bit_read.h"
#include "../bit_write.h"
#include "lr_string_mpz.h"
#include "string_utils.h"
#include "bit_read_string.h"
#define INF 0x7f7fffff

namespace Codecset
{

    uint8_t *encodeArray8_string(std::vector<std::string> &string_vec, int start_idx, const size_t length, uint8_t *res, size_t nvalue)
    {
        uint8_t *out = res;
        // vector of pointers to mpz_t
        std::vector<mpz_t *> ascii_vec;
        std::vector<int> index;
        for (int i = 0; i < length; i++)
        {
            auto *tmp_mpz = new mpz_t[1]; // TODO: no delete implemented; migrate to std::shared_ptr
            mpz_init(*tmp_mpz);
            ascii_vec.emplace_back(tmp_mpz);
        }
        for (int i = 0; i < length; i++)
        {
            // mpz_t tmp_val;
            // mpz_init(tmp_val);
            convertToASCII_mpz(string_vec[i + start_idx], ascii_vec[i]);

            // ascii_vec.emplace_back(&tmp_val);
            index.emplace_back(i);
        }

        string_mpz_lr mylr;
        mylr.caltheta(index, ascii_vec, length);
        // std::cout<<"theta0: "<<mylr.theta0<<" theta1: "<<mylr.theta1<<std::endl;
        std::vector<mpz_t *> delta;

        for (int i = 0; i < length; i++)
        {
            auto *tmp_mpz = new mpz_t[1]; // TODO: no delete implemented; migrate to std::shared_ptr
            mpz_init(*tmp_mpz);
            delta.emplace_back(tmp_mpz);
        }
        mpz_t max_delta;
        mpz_init(max_delta);
        for (auto i = 0; i < length; i++)
        {

            mpz_add(*delta[i], *delta[i], mylr.theta0);
            mpz_addmul_ui(*delta[i], mylr.theta1, i);

            mpz_sub(*delta[i], *ascii_vec[i], *delta[i]);

            if (mpz_cmpabs(*delta[i], max_delta) > 0)
            {
                mpz_set(max_delta, *delta[i]);
            }
        }

        uint32_t max_bit = 0;
        if (mpz_cmpabs_ui(max_delta, 0) > 0)
        {
            max_bit = bits_mpz(&max_delta) + 1;
        }
        mpz_clear(max_delta);

        memcpy(out, &max_bit, sizeof(uint32_t));
        out += sizeof(uint32_t);

        auto theta0_len = (mpz_sizeinbase(mylr.theta0, 2) + 7) / 8;
        memcpy(out, &theta0_len, sizeof(uint32_t));
        out += sizeof(uint32_t);
        mpz_export(out, NULL, -1, sizeof(uint8_t), 0, 0, mylr.theta0);
        out += theta0_len;

        auto theta1_len = (mpz_sizeinbase(mylr.theta1, 2) + 7) / 8;
        memcpy(out, &theta1_len, sizeof(uint32_t));
        out += sizeof(uint32_t);
        mpz_export(out, NULL, -1, sizeof(uint8_t), 0, 0, mylr.theta1);
        out += theta1_len;

        // for(auto item:delta){
        //     std::cout<<"delta: "<<*item<<std::endl;
        // }
        out = write_string_delta_string_mpz(delta, out, max_bit, length);
        mpz_clear(max_delta);
        for (int i = 0; i < length; i++)
        {
            mpz_clear(*ascii_vec[i]);
        }
        for (int i = 0; i < length; i++)
        {
            mpz_clear(*delta[i]);
        }

        return out;
    }

    void decodeArray8_string(uint8_t *in, const size_t length, size_t nvalue, std::vector<std::string> &string_vec)
    {

        uint32_t maxbits;
        uint8_t *tmpin = in;
        memcpy(&maxbits, tmpin, sizeof(uint32_t));
        // std::cout<<"max bit "<<maxbits<<std::endl;
        tmpin += sizeof(uint32_t);

        // double start = getNow();
        uint32_t theta0_len;
        memcpy(&theta0_len, tmpin, sizeof(uint32_t));
        tmpin += sizeof(uint32_t);
        // theta0 = reinterpret_cast<long_int*>(tmpin);
        mpz_t theta0;
        mpz_init(theta0);
        mpz_import(theta0, theta0_len, -1, 1, 0, 0, tmpin);
        tmpin += theta0_len;

        uint32_t theta1_len;
        memcpy(&theta1_len, tmpin, sizeof(uint32_t));
        tmpin += sizeof(uint32_t);
        mpz_t theta1;
        mpz_init(theta1);
        mpz_import(theta1, theta1_len, -1, 1, 0, 0, tmpin);
        tmpin += theta1_len;

        std::vector<mpz_t *> recover;
        for (int i = 0; i < length; i++)
        {
            auto *tmp_mpz = new mpz_t[1]; // TODO: no delete implemented; migrate to std::shared_ptr
            mpz_init(*tmp_mpz);
            recover.emplace_back(tmp_mpz);
        }

        read_all_bit_fix_string_mpz(tmpin, 0, 0, length, maxbits, &theta1, &theta0, recover);

        for (int i = 0; i < length; i++)
        {
            //std::cout << "out: " << *recover[i] << std::endl;
            std::string tmp_string = convertToString_mpz(recover[i]);
            //std::cout << "tmp_string: " << tmp_string << std::endl;
            string_vec.emplace_back(tmp_string);
        }
        // delete the mpz_t * array 
        for (int i = 0; i < length; i++)
        {
            mpz_clear(*recover[i]);
            delete recover[i];
        }
        // for (int i = 0; i < length; i++)
        // {
        //     mpz_clear(*recover[i]);

        // }
    }

} // namespace FastPFor

#endif /* SIMDFASTPFOR_H_ */
