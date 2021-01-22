/**
 * This code is released under the
 * Apache License Version 2.0 http://www.apache.org/licenses/.
 *
 * (c) Daniel Lemire, http://lemire.me/en/
 */
#ifndef PIECEWISE_H_
#define PIECEWISE_H_

#include "common.h"
#include "codecs.h"
#include "time.h"
#include "bit_opt.h"
#include "caltime.h"
#define INF 0x7f7fffff

namespace Codecset {

uint32_t lower_bound( uint32_t v,std::vector<uint32_t> segindex)
{
    uint32_t m;
    uint32_t x=0;
    uint32_t y=segindex.size()-1;
    while(x <= y)
    {
    
        m = x+(y-x)/2;
        if(v<segindex[m]) y = m-1;
        else x = m+1;
    }
    return y;
}

struct Segment{
    uint32_t start_index;
    float slope;
    uint8_t max_error; // how many bits we need to save this segment
    uint32_t start_byte;
    uint32_t start_key;
};
    
class piecewise : public IntegerCODEC {
public:
  using IntegerCODEC::encodeArray;
  using IntegerCODEC::decodeArray;
  using IntegerCODEC::randomdecodeArray;
  using IntegerCODEC::encodeArray8;
  using IntegerCODEC::decodeArray8;
  using IntegerCODEC::randomdecodeArray8;
   
  std::vector<uint32_t> segment_index;
  std::vector<Segment> q;
  int maxerror =63;

void setMaxError(int error){
    maxerror = error;
    std::cout<<"max error is set to: "<<maxerror<<std::endl;
}


float epsilon = 0.;
bool gt(float a, float b){
// greater than
    return (a > b + epsilon);
}
/*
uint8_t* write_delta(int *in,uint8_t* out, uint8_t l, int numbers){
    int code =0;
    int occupy = 0;

    uint8_t left_val = 0;

    for(int i = 0; i < numbers;i++){

        code = left_val;
        while(occupy<8){
            
            bool sign = 1;
            if (in[i] <= 0){
                sign = 0;
                in[i]=-in[i];
            }

            uint8_t value1= ((in[i] & ((1U<<(l-1))-1)) + (sign<<(l-1)));
            code += (value1<<occupy);
            occupy += l;

            i++;
            
        }//end while
        if(occupy>8){
            left_val = code >> 8;
            code = code & ((1U<<8) - 1);
        }
        else{
            left_val=0;
        
        }
        
        occupy = occupy % 8;
        
        out[0]= unsigned((uint8_t)code);
        out++;
        if(i==numbers){
            if (left_val){
                out[0]=(uint8_t)left_val;
                out++;
            }
            if(left_val==0&&occupy!=0){
                out[0]=(uint8_t)left_val;
                out++;
            }
            
        }
        i--;
    }

    
    return out;

}
*/
uint8_t* write_delta(int *in,uint8_t* out, uint8_t l, int numbers){
    int code =0;
    int occupy = 0;
    int endbit = (l*numbers);
    int end=0;
    uint8_t* start=out;
    if(endbit%8==0){
        end=endbit/8;
    }
    else{
        end = (int)endbit/8+1;
    }

    uint8_t* last=out+end;
    uint32_t left_val = 0;

    while(out<=last){

        
        while(occupy<8){
            
            bool sign = 1;
            if (in[0] <= 0){
                sign = 0;
                in[0]=-in[0];
            }
            
            uint32_t value1= ((in[0] & ((1U<<(l-1))-1)) + (sign<<(l-1)));
            //std::cout<<"add: "<<in[0]<<" value1 "<<value1<<std::endl;
            code += (value1<<occupy);
            occupy += l;

            in++;
            
        }//end while
        while(occupy>=8){
            left_val = code >> 8;
            //std::cout<<code<<std::endl;
            code = code & ((1U<<8) - 1);
            occupy-=8;
            out[0]= unsigned((uint8_t)code);
            code = left_val;
            //std::cout<<occupy<<" "<<left_val<<" "<<unsigned(out[0])<<std::endl;
            out++;
        }
    }
    

    
    return out;

}

uint8_t * encodeArray8(uint32_t *in, const size_t length,uint8_t *out, size_t nvalue) {
    uint8_t* start=out;
    std::vector<uint32_t> indexes;
    for(uint32_t i = 0; i < nvalue; i++){
        indexes.push_back(i);
    }   
    float high_slope = (float)INF;
    float low_slope = 0.;
    long long origin_key = in[0];
    int origin_index = indexes[0];
    int end_index = indexes[0];
    int start_byte = 0;
    int total_index =0;
    for (int i = 1; i < nvalue; i++){
        long long key = in[i];
        int id = indexes[i];
        float tmp_point_slope = ((key - origin_key)+0.0) / ((id - origin_index)+0.0);
        
        //if (tmp_point_slope >= low_slope && tmp_point_slope <= high_slope){
        if (gt(tmp_point_slope, low_slope) && gt(high_slope, tmp_point_slope)){
            float tmp_high_slope = ((key + maxerror - origin_key)+0.0) / ((id - origin_index)+0.0);
            float tmp_low_slope = ((key - maxerror - origin_key)+0.0) /((id - origin_index)+0.0);
            if (tmp_low_slope<0.0){
                //std::cout<<"zero!"<<std::endl;
                tmp_low_slope=0.0;
            }
            //std::cout<<tmp_high_slope<<","<<tmp_low_slope<<std::endl;
            
         //   if(tmp_high_slope<=high_slope){
           if(gt(high_slope, tmp_high_slope)){
                high_slope = tmp_high_slope;
            }
            //if(low_slope<=tmp_low_slope){
            if(gt(tmp_low_slope, low_slope)){
                low_slope = tmp_low_slope;
            }
            end_index = id;
            
        }
        else{
            
            float slope = (high_slope + low_slope) / 2.;
            int max_error = 0;

            if (end_index == origin_index){
                slope = 1.;
            }
            int seg_len = end_index - origin_index + 1;
            int * delta = static_cast<int *>(malloc(seg_len * sizeof(int)));
            
            
            for (int j = origin_index; j <= end_index; j++ ){
                long long  predict = (long long)in[origin_index] + (long long)(slope*(float)(j-origin_index));
                long long truth = (long long)in[j];
                int tmp = abs(predict-truth);
                
                delta[j-origin_index] = (int) truth-predict;
                total_index++;

                if (tmp > max_error){
                    max_error = tmp;
                    
                }
            }

            int tmp_bit = 0;
            if(max_error > 0.01){
                tmp_bit = ceil(log2(max_error))+2;
            }
            else{
                tmp_bit = 2;
            }
            Segment newseg;
            newseg.start_index = origin_index;
            newseg.slope = slope;
            newseg.max_error = tmp_bit;
            newseg.start_byte = start_byte;
            newseg.start_key = in[origin_index];
            //std::cout<<"bit_length: "<<tmp_bit<<" start: "<<origin_index<<" end: "<<end_index<<" slope: "<<slope<<std::endl;
            /*
            for(int k=0;k<(end_index - origin_index + 1);k++){
           std::cout<<k+origin_index<<" "<<delta[k]<<std::endl; 
            }
            */
            out = write_delta(delta, out, tmp_bit, (end_index - origin_index + 1));
            
            start_byte= out - start;
            
            q.push_back(newseg);
            segment_index.push_back(origin_index);
            free (delta);
            
            high_slope = (float)INF;
            low_slope = 0.0;
            origin_index = id;
            origin_key = key;
            end_index = id;   
            
        }

    }
    
    float slope = (high_slope + low_slope) / 2.;
    if (end_index == origin_index){
        slope = 1.;
    }
    
    int seg_len = end_index - origin_index + 1;
    int * delta = static_cast<int *>(malloc(seg_len * sizeof(int)));
    int max_error = 0;
    for (int j = origin_index; j <= end_index; j++ ){
        long long  predict = (long long)in[origin_index] + (long long)(slope*(float)(j-origin_index));
        long long truth = (long long)in[j];
        int tmp = abs(predict-truth);     
        delta[j-origin_index] = (int) truth-predict;

        total_index++;
        if (tmp > max_error){
            max_error = tmp;   
        }
    }

    uint8_t tmp_bit = 0;
    if(max_error > 0.01){
        tmp_bit = ceil(log2(max_error))+2;
    }
    else{
        tmp_bit = 2;
    }
    Segment newseg;
    newseg.start_index = origin_index;
    newseg.slope = slope;
    newseg.max_error = tmp_bit;
    newseg.start_byte = start_byte;
    newseg.start_key = in[origin_index];
    //std::cout<<"bit_length: "<<tmp_bit<<std::endl;
    out = write_delta(delta ,out, tmp_bit, (end_index - origin_index + 1));
    start_byte= out - start;
    q.push_back(newseg);
    segment_index.push_back(origin_index);
    free(delta);

    return out;
    
}
    
uint32_t *decodeArray8( uint8_t *in, const size_t length, uint32_t *out, size_t nvalue) {

    int start_byte=0;
    int end_index =0;
    int start_index = 0;
    int index=0;
    while(index<q.size()){
        Segment tmpseg = q[index];
        index++;
        start_byte = tmpseg.start_byte;
        start_index = tmpseg.start_index;
        if(index<q.size()){
            Segment tmpseg2 = q[index];
            end_index = tmpseg2.start_index-1;  
        }
        else{
            end_index = nvalue-1;
        }
        //std::cout<<"start byte: "<<start_byte<<" start_index: "<<start_index<<" end index: "<<end_index<<" start_key: "<<tmpseg.start_key<<std::endl;
        read_all_bit(in ,start_byte,start_index, (end_index-start_index+1),tmpseg.max_error,tmpseg.slope,tmpseg.start_key, out);
        
      
    }

    return out;
}
uint32_t randomdecodeArray8( uint8_t *in, const size_t l, uint32_t *out, size_t nvalue){
    
    //double start =getNow();
    Segment tmpseg = q[lower_bound(l,segment_index)];
    //double end = getNow();
    //double start2 = getNow();
    uint32_t tmp = read_bit(in ,tmpseg.max_error , l-tmpseg.start_index, tmpseg.slope,tmpseg.start_key, tmpseg.start_byte);
    //double end2 = getNow();
    //std::cout<<"find lower bound time: "<<(end-start)<<" read bit time: "<<(end2-start2)<<" rate is: "<<(end-start)/(end2-start2)<<std::endl;
    
    
    return tmp;

}

uint32_t* encodeArray( uint32_t *in, const size_t length, uint32_t *out,
                   size_t nvalue) {
    std::cout<<"Haven't implement. Please try uint8_t one..."<<std::endl;
    return out;
}
uint32_t *decodeArray( uint32_t *in, const size_t length,
                              uint32_t *out, size_t nvalue) {
    std::cout<<"Haven't implement. Please try uint8_t one..."<<std::endl;
    return out;
}
uint32_t randomdecodeArray(uint32_t *in, const size_t l,uint32_t *out, size_t nvalue){
    std::cout<<"Haven't implement. Please try uint8_t one..."<<std::endl;
    return 1;
}
uint32_t get_block_nums(){
     return q.size();
 }    
std::string name() const {
    return "FrameofReference"; 
}    
  
};



} // namespace FastPFor

#endif /* SIMDFASTPFOR_H_ */
