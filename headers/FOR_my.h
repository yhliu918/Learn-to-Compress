
#ifndef FOR_my_H_
#define FOR_my_H_

#include "common.h"
#include "codecs.h"
#include "bpacking.h"
#include "forutil.h"
#include "bit_write.h"
#include "bit_read.h"

namespace Codecset
{

    class FOR_my : public IntegerCODEC
    {
    public:
        using IntegerCODEC::decodeArray;
        using IntegerCODEC::decodeArray8;
        using IntegerCODEC::encodeArray;
        using IntegerCODEC::encodeArray8;
        using IntegerCODEC::init;
        using IntegerCODEC::randomdecodeArray;
        using IntegerCODEC::randomdecodeArray8;
        using IntegerCODEC::summation;

        int block_num;
        int block_size;

        void init(int blocks, int blocksize, int extra)
        {
            block_num = blocks;
            block_size = blocksize;
        }

        uint32_t *encodeArray(uint32_t *in, const size_t length, uint32_t *res, size_t nvalue)
        {

            return res;
        }
        uint32_t *decodeArray(uint32_t *in, const size_t length,
                              uint32_t *out, size_t nvalue)
        {
            uint32_t *res = out;
            nvalue = in[0];
            ++in;
            if (nvalue == 0)
                return in;
            uint32_t m = in[0];
            ++in;
            uint32_t M = in[0];
            ++in;

            int b = bits(static_cast<uint32_t>(M - m));
            if (b == 31)
            {
                b = 32;
            }
            // std::cout<<"bit "<<b<<" min "<<m<<" max "<<M<<std::endl;
#ifdef _OPENMP
#pragma omp parallel for
#endif
            if (b == 0)
            {
                for (int i = 0; i < (int)length; i++)
                {
                    out[i] = m;
                }
                return out;
            }
            for (uint32_t k = 0; k < length / 32; ++k)
            {
                unpack32[b](m, in + b * k, res + 32 * k);
            }
            res = res + length / 32 * 32;
            in = in + length / 32 * b;

            for (uint32_t k = length / 32 * 32; k + 16 <= length; k += 16, res += 16)
            {
                in = unpack16[b](m, in, res);
            }
            for (uint32_t k = length / 16 * 16; k + 8 <= length; k += 8, res += 8)
            {
                in = unpack8[b](m, in, res);
            }
            // we could pack the rest, but we don't  bother
            for (uint32_t k = length / 8 * 8; k < length; ++k, in++, res++)
            {
                res[0] = in[0];
            }
            return res;
        }
        uint32_t randomdecodeArray(uint32_t *in, const size_t l,
                                   uint32_t *out, size_t nvalue)
        {
            uint32_t tmpval = in[0];
            ++in;
            if (tmpval == 0)
                return in[0];
            uint32_t m = in[0];
            ++in;
            uint32_t M = in[0];
            ++in;
            int b = bits(static_cast<uint32_t>(M - m));
            if (b == 0)
            {
                return m;
            }
            if (b == 31)
            {
                b = 32;
            }
#ifdef _OPENMP
#pragma omp parallel for
#endif

            uint32_t recover = 0;
            in = in + ((int)l / 32) * b;
            int number_left = l - ((int)l / 32) * 32;
            // std::cout<<"nvalue "<<nvalue<<" number_left "<<number_left<<" nvalue - l "<<nvalue - l<<std::endl;
            // std::cout<<"to_find "<<l<<" block_size "<<block_size<<" bit "<<b<<std::endl;

            // if(((int)l/32) == ((int)nvalue/32)){

            //     /*
            //     if( number_left>=16  ){
            //         number_left -= 16;
            //         in = in +(int)ceil((double)b*16./32.);
            //         if(number_left>=8){
            //             number_left -= 8;
            //             in = in +(int)ceil((double)b*8./32.);
            //         }
            //         else{
            //             unpack8[b](m,in,out);
            //             return out[number_left];
            //         }

            //     }
            //     else if(number_left>=8){
            //         number_left -= 8;
            //         in = in +(int)ceil((double)b*8./32.);
            //         return in[number_left];

            //     }

            //     if(number_left>0){
            //         return in[number_left];
            //     }
            //     */
            //     uint32_t *res = new uint32_t[32];
            //     uint32_t *tmpres =res;
            //     for(uint32_t k=nvalue/32*32; k+16<=nvalue; k+=16,res+=16) {
            //         in = unpack16[b](m,in,res);
            //     }
            //     for(uint32_t k=nvalue/16*16; k+8<=nvalue; k+=8,res+=8) {
            //         in = unpack8[b](m,in,res);
            //     }
            // // we could pack the rest, but we don't  bother
            //     for(uint32_t k=nvalue/8*8; k<nvalue; ++k,in++,res++) {
            //         res[0] = in [0];
            //     }
            //     recover = tmpres[number_left];
            //     free(tmpres);
            //     return recover;
            // }
            // else{
            uint32_t number_occupy = (number_left * b) / 32;
            in += number_occupy;

            if (b == 32)
            {
                return in[0] + m;
            }

            long long bit_left = number_left * b - number_occupy * 32;

            if (32 - bit_left >= b)
            {
                recover = (in[0] >> bit_left) & ((1U << b) - 1);
                recover += m;

                return recover;
            }
            else
            {
                recover = ((in[1] & (((1U << (b + bit_left - 32)) - 1))) << (32 - bit_left)) + (in[0] >> bit_left);
                recover += m;
                return recover;
            }
        }

        uint8_t *encodeArray8(uint32_t *in, const size_t length, uint8_t *res,
                              size_t nvalue)
        {
            uint8_t *out = res;
            uint32_t max_record = 0;
            uint32_t min_record = UINT32_MAX;
            for (int i = 0; i < length; i++)
            {
                if (in[i] > max_record)
                    max_record = in[i];
                if (in[i] < min_record)
                    min_record = in[i];
            }
            uint32_t max_delta = max_record - min_record;
            int tmp_bits = 0;
            if (max_delta)
            {
                tmp_bits = bits(max_delta);
            }
            std::vector<int> delta;
            for (int i = 0; i < length; i++)
            {
                delta.push_back(in[i] - min_record);
            }
            memcpy(out, &tmp_bits, 1);
            out += 1;
            if (tmp_bits >= sizeof(uint32_t) * 8 - 1)
            {
                for (auto i = 0; i < length; i++)
                {
                    memcpy(out, &in[i], sizeof(uint32_t));
                    out += sizeof(uint32_t);
                }
                return out;
            }

            memcpy(out, &min_record, sizeof(uint32_t));
            out += sizeof(uint32_t);
            memcpy(out, &max_record, sizeof(uint32_t));
            out += sizeof(uint32_t);

            if (tmp_bits)
            {

                out = write_FOR_int_T(delta.data(), out, tmp_bits, length);
            }
            return out;
        }

        uint32_t *decodeArray8(uint8_t *in, const size_t length,
                               uint32_t *out, size_t nvalue)
        {

            uint8_t maxerror;
            uint8_t *tmpin = in;
            memcpy(&maxerror, tmpin, 1);
            tmpin++;
            if (maxerror >= sizeof(uint32_t) * 8 - 1)
            {
                memcpy(out, tmpin, length * sizeof(uint32_t));
                return out;
            }
            uint32_t base = 0;
            memcpy(&base, tmpin, sizeof(base));
            tmpin += sizeof(base);
            uint32_t max = 0;
            memcpy(&max, tmpin, sizeof(max));
            tmpin += sizeof(max);

            if (maxerror)
            {

                read_all_bit_FOR<uint32_t>(tmpin, 0, length, maxerror, base, out);
            }
            else
            {
                for (int i = 0; i < length; i++)
                {
                    out[i] = base;
                }
            }

            return out;
        }

        uint32_t randomdecodeArray8(uint8_t *in, const size_t l, uint32_t *out, size_t nvalue)
        {
            uint8_t *tmpin = in;
            uint8_t maxbits;
            memcpy(&maxbits, tmpin, sizeof(uint8_t));
            tmpin += sizeof(uint8_t);
            if (maxbits >= sizeof(uint32_t) * 8 - 1)
            {
                uint32_t tmp_val = reinterpret_cast<uint32_t *>(tmpin)[l];
                return tmp_val;
            }

            uint32_t m;
            memcpy(&m, tmpin, sizeof(m));
            tmpin += sizeof(m);
            uint32_t max;
            memcpy(&max, tmpin, sizeof(max));
            tmpin += sizeof(max);

            uint32_t tmp_val = m;
            if (maxbits)
            {
                tmp_val += read_FOR_int<uint32_t>(tmpin, maxbits, l);
            }
            return tmp_val;
        }
        uint64_t summation(uint8_t *in, const size_t l, size_t nvalue)
        {
            uint32_t *res = (uint32_t *)malloc(nvalue * sizeof(uint32_t));
            decodeArray8(in, nvalue, res, nvalue);
            uint64_t sum = 0;
            for (int i = 0; i < (int)nvalue; i++)
            {
                sum += res[i];
            }
            return sum;
        }
        uint32_t get_block_nums()
        {
            return 1;
        }
        void destroy() {}
        std::string name() const
        {
            return "FOR_my";
        }
    };

} // namespace FastPFor

#endif /* SIMDFASTPFOR_H_ */
