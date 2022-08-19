
#ifndef delta_my_H_
#define delta_my_H_

#include "common.h"
#include "codecs.h"

namespace Codecset
{

    class delta_my : public IntegerCODEC
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
        // max_bit+number1+delta
        uint32_t *encodeArray(uint32_t *in, const size_t length, uint32_t *res, size_t nvalue)
        {
            uint8_t * out=reinterpret_cast<uint8_t*>(res);
            uint8_t * mark_out = out;
            out = encodeArray8(in,length,out,nvalue);
            res = reinterpret_cast<uint32_t*>(mark_out);
            uint32_t *tmp_res = reinterpret_cast<uint32_t*>(out);
            return tmp_res;
        }

        uint32_t *decodeArray(uint32_t *in, const size_t length, uint32_t *out, size_t nvalue)
        {
            uint8_t *tmpin = reinterpret_cast<uint8_t *>(in);
            return decodeArray8(tmpin, length, out, nvalue);
        }

        uint32_t randomdecodeArray(uint32_t *in, const size_t l, uint32_t *out, size_t nvalue)
        {
            std::cout<<"This method doesn't support random access."<<std::endl;
            return 0;
        }

        uint8_t *encodeArray8(uint32_t *in, const size_t length, uint8_t *res, size_t nvalue)
        {
            uint8_t *out = res;
            int *delta = new int[length];
            int max_delta=0;
            for(uint32_t i=0;i<length-1;i++){
                delta[i] = in[i+1]-in[i];
                if(abs(delta[i])>max_delta){
                    max_delta =  abs(delta[i]);
                }
            }
            int tmp_bit = bits(max_delta)+1;
            // std::cout<<tmp_bit<<std::endl;
            if(max_delta==0){
                tmp_bit=0;
            }
            
            out[0]=(uint8_t)tmp_bit;
            out++;
            out = write_delta_default(&in[0],out, 32, 1);
            if(tmp_bit==0){
                free(delta);
                return out;
            }
            if(tmp_bit>=32){
                out = write_delta_default(in,out,32,length);
            }
            else{
                out = write_delta(delta, out, tmp_bit, length-1);
            }
    
            free(delta);
            return out;
        }

        uint32_t *decodeArray8(uint8_t *in, const size_t length, uint32_t *out, size_t nvalue)
        {
            uint32_t* res = out;
            //start_index + bit + theta0 + theta1 + numbers + delta
            uint8_t* tmpin = in;
            uint8_t maxerror = tmpin[0];
            tmpin++;
            uint32_t base;
            memcpy(&base, tmpin, sizeof(base));
            tmpin += sizeof(uint32_t);
            if (maxerror) {
                res[0] = base;
                read_all_bit_Delta<uint32_t>(tmpin, 0, length - 1, maxerror, base, res + 1);
            }
            else {
                for (int j = 0;j < length;j++) {
                    res[j] = base;
                }
            }

            return out;

        }
        uint32_t randomdecodeArray8(uint8_t *in, const size_t l, uint32_t *out, size_t nvalue)
        {
            uint8_t* tmpin = in;
            uint8_t maxbits;
            memcpy(&maxbits, tmpin, sizeof(uint8_t));
            tmpin += sizeof(uint8_t);

            if (maxbits == sizeof(uint32_t) * 8) {
                uint32_t tmp_val = reinterpret_cast<uint32_t*>(tmpin)[l];
                return tmp_val;
            }

            uint32_t base;
            memcpy(&base, tmpin, sizeof(base));
            tmpin += sizeof(base);

            uint32_t tmp_val = base;
            if (maxbits) {
                tmp_val = read_Delta_int(tmpin, maxbits, l, base);
            }
            return tmp_val;

        }
        uint64_t summation( uint8_t *in, const size_t l, size_t nvalue){
    
            return 0;
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
            return "delta_my";
        }
    };
}

#endif