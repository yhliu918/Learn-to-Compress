
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
        uint8_t *encodeArray8_int(T *data, const size_t length, uint8_t *res, size_t nvalue)
        {
            uint8_t *out = res;

            std::vector<T> delta;
            std::vector<bool> signvec;
            T max_error = 0;

            for (auto i = 0; i < length-1; i++)
            {
                T tmp_val;
                if ( data[i+1] > data[i])
                {
                    tmp_val = data[i+1] - data[i];
                    signvec.emplace_back(true); // means positive
                }
                else
                {
                    tmp_val = data[i] -  data[i+1];
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
            
            if(max_bit>sizeof(T)*8){
                max_bit = sizeof(T)*8;
            }
            memcpy(out, &max_bit, sizeof(max_bit));
            out += sizeof(max_bit);
            if(max_bit== sizeof(T)*8){
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
                out = write_delta_int_T(delta.data(),signvec, out, max_bit, length);
            }

            return out;
        }


        T randomdecodeArray8(uint8_t *in, int to_find, uint32_t *out, size_t nvalue)
        {
            
            uint8_t *tmpin = in;
            uint8_t maxbits;
            memcpy(&maxbits, tmpin, sizeof(uint8_t));
            tmpin += sizeof(uint8_t);

            if(maxbits==sizeof(T)*8){
                T tmp_val = reinterpret_cast<T *>(tmpin)[to_find];
                return tmp_val;
            }

            T base;
            memcpy(&base, tmpin, sizeof(T));
            tmpin += sizeof(T);
            
            T tmp_val = base;
            if(maxbits){
                tmp_val = read_Delta_int(tmpin, maxbits, to_find, base);
            }
            return tmp_val;
        }



    };

} // namespace FastPFor

#endif /* SIMDFASTPFOR_H_ */