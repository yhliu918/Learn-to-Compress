
#ifndef delta_H_
#define delta_H_

#include "common.h"
#include "codecs.h"

namespace Codecset
{

    class delta : public IntegerCODEC
    {
    public:
        using IntegerCODEC::decodeArray;
        using IntegerCODEC::decodeArray8;
        using IntegerCODEC::encodeArray;
        using IntegerCODEC::encodeArray8;
        using IntegerCODEC::init;
        using IntegerCODEC::randomdecodeArray;
        using IntegerCODEC::randomdecodeArray8;

        int block_num;
        int block_size;

        void init(int blocks, int blocksize, int extra)
        {
            block_num = blocks;
            block_size = blocksize;
        }

        uint32_t *encodeArray(uint32_t *in, const size_t length, uint32_t *res, size_t nvalue)
        {
            uint32_t *out = res;
            uint32_t M = 0;
            for (uint32_t i = 1; i < length; ++i)
            {
                if (in[i] - in[i - 1] > M)
                    M = in[i] - in[i - 1];
            }
            *out = in[0];
            ++out;
            int b = bits(M);

            *out = static_cast<uint32_t>(b);
            ++out;

            if (b == 0)
                return out;

            int lbit = 32; //当前32位数还剩下的位数
            *out = 0;
            for (uint32_t i = 1; i < length; ++i)
            {
                if (b <= lbit)
                {
                    *out = (*out) << b;
                    *out |= in[i] - in[i - 1];
                    lbit -= b;
                }
                else
                {
                    *out = (*out) << lbit;
                    uint32_t value = (in[i] - in[i - 1]) >> (b - lbit); //只保留高lbit位
                    *out |= value;
                    ++out;    //当前32位数已经用满，用下一个32位
                    *out = 0; //initialize
                    int mask = (1 << (b - lbit)) - 1;
                    *out |= (in[i] - in[i - 1]) & mask;
                    lbit = 32 - (b - lbit);
                }
            }
            if (lbit < 32)
            {
                *out = (*out) << lbit;
                ++out;
                *out = 0;
            }
            return out;
        }

        uint32_t *decodeArray(uint32_t *in, const size_t length, uint32_t *out, size_t nvalue)
        {
            out[0] = *in;
            ++in;
            int b = static_cast<int>(*in);
            ++in;
            if (b == 0)
            {
                for (uint32_t i = 1; i < length; ++i)
                {
                    out[i] = out[0];
                }
                return out + length - 1;
            }
            int lbit = 32;
            int mask = (1 << b) - 1;
            for (uint32_t i = 1; i < length; ++i)
            {
                if (b <= lbit)
                {
                    out[i] = (*in) >> (lbit - b);
                    out[i] = out[i] & mask;
                    out[i] += out[i - 1];
                    lbit -= b;
                }
                else
                {
                    int tmp_mask = (1 << lbit) - 1;
                    out[i] = (*in) & tmp_mask;
                    out[i] = out[i] << (b - lbit);
                    ++in;
                    out[i] |= (*in) >> (32 - (b - lbit));
                    out[i] += out[i - 1];
                    lbit = 32 - (b - lbit);
                }
            }
            return out + length - 1;
        }

        uint32_t randomdecodeArray(uint32_t *in, const size_t l, uint32_t *out, size_t nvalue)
        {
            uint32_t M = *in;
            ++in;
            int b = static_cast<int>(*in);
            ++in;
            if (b == 0 || l == 0)
                return M;
            uint32_t recover = M;
            for (uint32_t i = 1; i <= l; ++i)
            {
                int dis = (b * (i-1)) / 32;
                uint32_t *tmp_des = in + dis;
                int st = 32 - (b * (i-1)) % 32;

                uint32_t tmp_res = 0;
                if (st < b)
                {
                    int mask = (1 << st) - 1;
                    tmp_res = (*tmp_des) & mask;
                    tmp_res = tmp_res << (b - st);
                    tmp_des = tmp_des + 1;
                    mask = (1 << (b - st)) - 1;
                    uint32_t tmp_val = (*tmp_des) >> (32 - b + st);
                    tmp_val &= mask;
                    tmp_res |= tmp_val;
                }
                else
                {
                    tmp_res = (*tmp_des) >> (st - b);
                    int mask = (1 << b) - 1;
                    tmp_res =tmp_res & mask;
                }
                recover += tmp_res;
            }
            return recover;
        }

        uint8_t *encodeArray8(uint32_t *in, const size_t length, uint8_t *res, size_t nvalue)
        {
            uint32_t *out = reinterpret_cast<uint32_t *>(res);
            uint32_t *mark_out = out;
            out = encodeArray(in, length, out, nvalue);

            res = reinterpret_cast<uint8_t *>(mark_out);
            uint8_t *tmp_des = reinterpret_cast<uint8_t *>(out);
            return tmp_des;
        }

        uint32_t *decodeArray8(uint8_t *in, const size_t length, uint32_t *out, size_t nvalue)
        {
            uint32_t *tmpin = reinterpret_cast<uint32_t *>(in);
            return decodeArray(tmpin, length, out, nvalue);
        }
        uint32_t randomdecodeArray8(uint8_t *in, const size_t l, uint32_t *out, size_t nvalue)
        {
            uint32_t *tmpin = reinterpret_cast<uint32_t *>(in);
            uint32_t tmp = randomdecodeArray(tmpin, l, out, nvalue);
            return tmp;
        }
        uint32_t get_block_nums()
        {
            return 1;
        }
        void destroy()
        {
        }
        std::string name() const
        {
            return "delta";
        }
    };
}

#endif