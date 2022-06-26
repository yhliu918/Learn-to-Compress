
#ifndef FOR_INTEGER_TEMPLATE_H_
#define FOR_INTEGER_TEMPLATE_H_

#include "common.h"
#include "codecs.h"
#include "bit_write.h"
#include "bit_read.h"
#define INF 0x7f7fffff

namespace Codecset
{
    template <typename T>
    class FOR_int
    {
    public:
        uint8_t *encodeArray8_int(T *data, const size_t length, uint8_t *res, size_t nvalue)
        {
            uint8_t *out = res;

            
            T m = data[0];
            T M = data[0];
            for(uint32_t i = 1; i < length; ++i) {
                if(data[i]>M) M=data[i];
                if(data[i]<m) m=data[i];
            }

            std::vector<T> delta;

            for (auto i = 0; i < length; i++)
            {
                T tmp_val = data[i] - m;
                delta.emplace_back(tmp_val);
            }


            T max_error = M - m;
            uint8_t max_bit = 0;
            
            if (max_error)
            {
                max_bit = bits_int_T(max_error); // without sign bit
            }
            if(max_bit>sizeof(T)*8){
                max_bit = sizeof(T)*8;
            }
            

            memcpy(out, &max_bit, sizeof(uint8_t));
            out += sizeof(uint8_t);

            if(max_bit== sizeof(T)*8){
                for (auto i = 0; i < length; i++)
                {
                    memcpy(out, &data[i], sizeof(T));
                    out += sizeof(T);
                }
                return out;
            }

            memcpy(out, &m, sizeof(m));
            out += sizeof(m);
            memcpy(out, &M, sizeof(M));
            out += sizeof(M);
            if (max_bit)
            {
                out = write_FOR_int_T(delta.data(), out, max_bit, length);
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


            T m;
            memcpy(&m, tmpin, sizeof(m));
            tmpin += sizeof(m);

            T M;
            memcpy(&M, tmpin, sizeof(M));
            tmpin += sizeof(M);
            
            
            T tmp_val = m;
            if(maxbits){
                tmp_val += read_FOR_int<T>(tmpin, maxbits, to_find);
            }
            return tmp_val;
        }



    };

} // namespace FastPFor

#endif /* SIMDFASTPFOR_H_ */