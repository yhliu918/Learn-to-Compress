
#ifndef DELTA_INTEGER_TEMPLATE_H_
#define DELTA_INTEGER_TEMPLATE_H_

#include "common.h"
#include "codecs.h"
#include "bit_write.h"
#include "lr.h"
#include "bit_read.h"
#define INF 0x7f7fffff

namespace Codecset
{
    template <typename T>
    class Delta_int
    {
    public:
        std::vector<std::pair<int,int>> mul_add_diff_set;
        int blocks;
        int block_size;
        
        void init(int block, int block_s){
            blocks = block;
            block_size = block_s;
        }
        uint8_t* encodeArray8_int(T* data, const size_t length, uint8_t* res, size_t nvalue)
        {
            uint8_t* out = res;

            std::vector<T> delta;
            std::vector<bool> signvec;
            T max_error = 0;

            for (auto i = 0; i < length - 1; i++)
            {
                T tmp_val;
                if (data[i + 1] > data[i])
                {
                    tmp_val = data[i + 1] - data[i];
                    signvec.emplace_back(true); // means positive
                }
                else
                {
                    tmp_val = data[i] - data[i + 1];
                    signvec.emplace_back(false); // means negative
                }

                delta.emplace_back(tmp_val);

                if (tmp_val > max_error)
                {
                    max_error = tmp_val;
                }
            }



            uint8_t max_bit = 0;
            if (max_error)
            {
                max_bit = bits_int_T(max_error) + 1;
            }

            if (max_bit > sizeof(T) * 8) {
                max_bit = sizeof(T) * 8;
            }
            memcpy(out, &max_bit, sizeof(max_bit));
            out += sizeof(max_bit);
            if (max_bit == sizeof(T) * 8) {
                for (auto i = 0; i < length; i++)
                {
                    memcpy(out, &data[i], sizeof(T));
                    out += sizeof(T);
                }
                return out;
            }

            memcpy(out, &data[0], sizeof(T));
            out += sizeof(T);

            if (max_bit)
            {
                out = write_delta_int_T(delta, signvec, out, max_bit, length);
            }

            return out;
        }

        T* decodeArray8(uint8_t* in, const size_t length, T* out, size_t nvalue) {
            T* res = out;
            //start_index + bit + theta0 + theta1 + numbers + delta
            uint8_t* tmpin = in;
            uint8_t maxerror = tmpin[0];
            tmpin++;
            if (maxerror == 127) {
                T tmp_val;
                memcpy(&tmp_val, tmpin, sizeof(tmp_val));
                res[0] = tmp_val;
                res++;
                return out;
            }
            if (maxerror == 126) {
                T tmp_val;
                memcpy(&tmp_val, tmpin, sizeof(tmp_val));
                res[0] = tmp_val;
                res++;
                memcpy(&tmp_val, tmpin + sizeof(T), sizeof(tmp_val));
                res[0] = tmp_val;
                res++;
                return out;
            }

            T base;
            memcpy(&base, tmpin, sizeof(T));
            tmpin += sizeof(T);
            if (maxerror) {
                res[0] = base;
                read_all_bit_Delta<T>(tmpin, 0, length - 1, maxerror, base, res + 1);
            }
            else {
                for (int j = 0;j < length;j++) {
                    res[j] = base;
                }
            }

            return out;
        }




        T randomdecodeArray8(uint8_t* in, int to_find, uint32_t* out, size_t nvalue)
        {

            uint8_t* tmpin = in;
            uint8_t maxbits;
            memcpy(&maxbits, tmpin, sizeof(uint8_t));
            tmpin += sizeof(uint8_t);

            if (maxbits == sizeof(T) * 8) {
                T tmp_val = reinterpret_cast<T*>(tmpin)[to_find];
                return tmp_val;
            }

            T base;
            memcpy(&base, tmpin, sizeof(T));
            tmpin += sizeof(T);

            T tmp_val = base;
            if (maxbits) {
                tmp_val = read_Delta_int(tmpin, maxbits, to_find, base);
            }
            return tmp_val;
        }



    };

} // namespace FastPFor

#endif /* SIMDFASTPFOR_H_ */
