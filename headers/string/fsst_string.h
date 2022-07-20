
#ifndef FSST_STRING_TEMPLATE_H_
#define FSST_STRING_TEMPLATE_H_

#include "../common.h"
#include "../codecs.h"
#include "fsst.h"

#define INF 0x7f7fffff
#include <limits>
using namespace std;
namespace Codecset
{

    class FSST_string {
    public:
        uint64_t encodeArray8_string(const std::vector<std::string>& data, std::vector<uint32_t>& offset) {
            compressedData.clear();
            offsets.clear();

            vector<unsigned long> rowLens, compressedRowLens;
            vector<unsigned char*> rowPtrs, compressedRowPtrs;
            rowLens.reserve(data.size());
            compressedRowLens.resize(data.size());
            rowPtrs.reserve(data.size());
            compressedRowPtrs.resize(data.size() + 1);
            unsigned totalLen = 0;
            for (auto& d : data) {
                totalLen += d.size();
                rowLens.push_back(d.size());
                rowPtrs.push_back(reinterpret_cast<unsigned char*>(const_cast<char*>(d.data())));
            }

            auto startTime = std::chrono::steady_clock::now();
            auto encoder = fsst_create(data.size(), rowLens.data(), rowPtrs.data(), false);
            auto createTime = std::chrono::steady_clock::now();
            vector<unsigned char> compressionBuffer;
            compressionBuffer.resize(16 + 2 * totalLen);
            auto compressTime = std::chrono::steady_clock::now();
            fsst_compress(encoder, data.size(), rowLens.data(), rowPtrs.data(), compressionBuffer.size(), compressionBuffer.data(), compressedRowLens.data(), compressedRowPtrs.data());
            auto stopTime = std::chrono::steady_clock::now();
            unsigned long compressedLen = data.empty() ? 0 : (compressedRowPtrs[data.size() - 1] + compressedRowLens[data.size() - 1] - compressionBuffer.data());

            compressedData.resize(compressedLen + 8192);
            memcpy(compressedData.data(), compressionBuffer.data(), compressedLen);
            offsets.reserve(data.size());
            compressedRowPtrs[data.size()] = compressionBuffer.data() + compressedLen;
            for (unsigned index = 0, limit = data.size(); index != limit; ++index)
                offsets.push_back(compressedRowPtrs[index + 1] - compressionBuffer.data());
            uint64_t result = compressedData.size()/* + (offsets.size() * sizeof(unsigned))*/;
            {
                unsigned char buffer[sizeof(fsst_decoder_t)];
                unsigned dictLen = fsst_export(encoder, buffer);
                fsst_destroy(encoder);
                result += dictLen;

                fsst_import(&decoder, buffer);
            }
            offset = offsets;
            offsets.clear();
            return result;
        }

        /// Decompress a single row. The target buffer is guaranteed to be large enough


        std::string randomdecode_string(int idx, int offset_val, int offset_val2) {
            vector<char> target;
            target.resize((offset_val2-offset_val)*10);
            char* writer = target.data();
            auto limit = writer + target.size();
            auto data = compressedData.data();
            auto start = offset_val;
            auto end = offset_val2;
            unsigned len = fsst_decompress(&decoder, end - start, data + start, limit - writer, reinterpret_cast<unsigned char*>(writer));
            return  std::string(writer, writer + len);
        }


    private:
        fsst_decoder_t decoder;
        /// The compressed data
        vector<unsigned char> compressedData;
        /// The offsets
        vector<unsigned> offsets;
        uint8_t* compressedOffset;
    };

} // namespace FastPFor

#endif /* SIMDFASTPFOR_H_ */
