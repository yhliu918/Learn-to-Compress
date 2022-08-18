/**
 * This code is released under the
 * Apache License Version 2.0 http://www.apache.org/licenses/.
 *
 * (c) Daniel Lemire, http://lemire.me/en/
 */
#ifndef MASKVBYTE_H_
#define MASKVBYTE_H_

#include "common.h"
#include "codecs.h"
#include "varintencode.h"
#include "varintdecode.h"


namespace Codecset {


    class MaskVByte : public IntegerCODEC {
    public:
        using IntegerCODEC::encodeArray;
        using IntegerCODEC::decodeArray;
        using IntegerCODEC::randomdecodeArray;
        using IntegerCODEC::encodeArray8;
        using IntegerCODEC::decodeArray8;
        using IntegerCODEC::randomdecodeArray8;
        using IntegerCODEC::init;
        using IntegerCODEC::summation;

        int block_num;
        int block_size;

        void init(int blocks, int blocksize, int extra) {
            block_num = blocks;
            block_size = blocksize;
        }


        uint32_t* encodeArray(uint32_t* in, const size_t length, uint32_t* out, size_t nvalue) {
            std::cout << "Haven't implement. Please try encodeArray8 ..." << std::endl;
            return out;

        }
        uint32_t* decodeArray(uint32_t* in, const size_t length,
            uint32_t* out, size_t nvalue) {
            std::cout << "Haven't implement. Please try decodeArray8 ..." << std::endl;
            return out;

        }
        uint32_t randomdecodeArray(uint32_t* in, const size_t l,
            uint32_t* out, size_t nvalue) {
            std::cout << "Haven't implement. Cannot random access..." << std::endl;
            return 1;

        }

        uint8_t* encodeArray8(uint32_t* in, const size_t length, uint8_t* out,
            size_t nvalue) {

            return out + vbyte_encode(in, length, out);

        }
        uint32_t* decodeArray8(uint8_t* in, const size_t length,
            uint32_t* out, size_t nvalue) {
            masked_vbyte_decode(in, out, length);
            return out;
        }
        uint32_t randomdecodeArray8(uint8_t* in, const size_t l, uint32_t* out, size_t nvalue) {
            std::cout << "Haven't implement. Cannot random access...." << std::endl;
            return 1;
        }
        uint64_t summation(uint8_t* in, const size_t l, size_t nvalue) {
            return 0;
        }
        uint32_t get_block_nums() {
            return 1;
        }
        std::string name() const {
            return "MaskVByte";
        }
        void destroy() {}
    };



} // namespace FastPFor

#endif /* SIMDFASTPFOR_H_ */
