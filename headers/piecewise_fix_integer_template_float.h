
#ifndef LECO_FLOAT_INTEGER_TEMPLATE_H_
#define LECO_FLOAT_INTEGER_TEMPLATE_H_

#include "common.h"
#include "codecs.h"
#include "bit_write.h"
#include "lr.h"
#include "bit_read.h"
#define INF 0x7f7fffff

namespace Codecset
{
    template <typename T>
    class Leco_int_float
    {
    public:
        uint8_t *encodeArray8_int(T *data, const size_t length, uint8_t *res, size_t nvalue)
        {
            uint8_t *out = res;
            std::vector<double> indexes;
            std::vector<double> keys;
            for(uint32_t i = 0; i < length; i++){
                indexes.emplace_back((double) i);
                keys.emplace_back((double) data[i]);
            }

            lr mylr;
            mylr.caltheta(indexes,keys,length);
            T theta0 = mylr.theta0;
            float theta1 = (float)mylr.theta1;

            std::vector<T> delta;
            std::vector<bool> signvec;
            T max_error =0;

            for (auto i = 0; i < length; i++)
            {
                T tmp_val;
                if ( data[i] > (T)((double)theta0+theta1*(float)i))
                {
                    tmp_val = data[i] - (T)((double)theta0+theta1*(float)i);
                    signvec.emplace_back(true); // means positive
                }
                else
                {
                    tmp_val =(T)((double)theta0+theta1*(float)i) -  data[i];
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

            memcpy(out, &theta0, sizeof(theta0));
            out += sizeof(theta0);
            memcpy(out, &theta1, sizeof(float));
            out += sizeof(float);
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

            T theta0;
            memcpy(&theta0, tmpin, sizeof(theta0));
            tmpin += sizeof(theta0);

            float theta1;
            memcpy(&theta1, tmpin, sizeof(float));
            tmpin += sizeof(float);
            
            
            T tmp_val;
            if(maxbits){
                tmp_val = read_bit_fix_int_float<T>(tmpin, maxbits, to_find, theta1, theta0);
            }
            else{
                tmp_val = ((double)theta0 + (float)to_find * theta1);
            }
            return tmp_val;
        }



    };

} // namespace FastPFor

#endif /* SIMDFASTPFOR_H_ */
