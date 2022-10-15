
#include "codecs.h"
#include "codecfactory.h"
#include "FOR.h"
#include "rle.h"
#include "piecewise_fix.h"
#include "nonlinear_fix.h"

namespace Codecset
{
    template <class Codec1, class Codec2>
    class CombinedCodec : public IntegerCODEC
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
        
        CombinedCodec() : codec1(), codec2() {}
        Codec1 codec1;
        Codec2 codec2;

        void init(int blocks, int blocksize, int extra)
        {
            block_num = blocks;
            block_size = blocksize;
        }

        
        uint8_t *encodeArray8(uint32_t *in, const size_t length,uint8_t *out, size_t nvalue) 
        {

            uint8_t * tmp_out = out;
            uint32_t compress_f1_size = 0;
            uint8_t *descriptor_tmp = (uint8_t *)malloc(block_size * sizeof(uint64_t) * 2);
            //memset(descriptor_tmp,0,sizeof(descriptor_tmp));
            uint8_t *res = descriptor_tmp;
            res = codec1.encodeArray8(in , block_size, descriptor_tmp, nvalue);
            compress_f1_size = res-descriptor_tmp;
            int resize_size =(int) ceil((double)compress_f1_size/4.0);
            descriptor_tmp=(uint8_t*)realloc(descriptor_tmp, resize_size*4);
            uint32_t * descriptor_tmp32=reinterpret_cast<uint32_t*>(descriptor_tmp);
            tmp_out = write_delta_default(&compress_f1_size,out, 32, 1);
            tmp_out = codec2.encodeArray8(descriptor_tmp32, resize_size, tmp_out, nvalue);
            free(descriptor_tmp);

            /*
            int block_num2 = ceil((double)resize_size/(double)block_size2);
            int bias = 0;
            for(int i=0;i<block_num2-1;i++){
                tmp_out = f2->encodeArray8(descriptor_tmp+i*block_size2, block_size2, out+bias, nvalue);
                bias = tmp_out - out;
            }
            tmp_out = f2->encodeArray8(descriptor_tmp+(block_num2-1)*block_size2, (resize_size - (block_num2-1)*block_size2), out+bias, nvalue);
            */

            return tmp_out;                         
                                   
             
        }
        


        uint32_t *decodeArray8(uint8_t *in, const size_t length,uint32_t *out, size_t nvalue)
        {
            uint32_t code = in[0];
            code += (in[1]<<8);
            code += (in[2]<<16);
            code += (in[3]<<24);
            in+=4;
            int resize_size =(int) ceil((double)code/4.0);
            uint32_t * recover1 = (uint32_t*)malloc(resize_size * sizeof(uint64_t));
            codec2.decodeArray8(in, resize_size, recover1, nvalue); 
            uint8_t * recover_tmp=reinterpret_cast<uint8_t*>(recover1);
            recover_tmp=(uint8_t*)realloc(recover_tmp, code);
            codec1.decodeArray8(recover_tmp, block_size, out, nvalue);
            free(recover1);
            
            return out; 
        }
        uint32_t randomdecodeArray8(uint8_t *in, const size_t l,uint32_t *out, size_t nvalue) {
            uint32_t code = in[0];
            code += (in[1]<<8);
            code += (in[2]<<16);
            code += (in[3]<<24);
            in+=4;
            int resize_size =(int) ceil((double)code/4.0);
            uint32_t * recover1 = (uint32_t*)malloc(resize_size * sizeof(uint64_t));
            codec2.decodeArray8(in, resize_size, recover1, nvalue); 
            uint8_t * recover_tmp=reinterpret_cast<uint8_t*>(recover1);
            uint32_t tmp_val = codec1.randomdecodeArray8(recover_tmp, l, out, nvalue);
            free(recover1);
            
            
            return tmp_val; 
        }

        uint32_t *decodeArray(uint32_t *in, const size_t length,
                              uint32_t *out, size_t nvalue) { return 0; }
        uint32_t randomdecodeArray(uint32_t *in, const size_t l,
                                   uint32_t *out, size_t nvalue) { return 0; }
        uint32_t *encodeArray(uint32_t *in, const size_t length,uint32_t *out, size_t nvalue) 
        {                      
            std::cout<<"Haven't implement. Please try uint8_t one..."<<std::endl;                 
            return 0; 
        }
        uint64_t summation( uint8_t *in, const size_t l, size_t nvalue){
            return 0;
        }
        uint32_t get_block_nums() { return block_num; }

        void destroy() {}

        std::string name() const { return "Combined"; }

        
        
    };
}
