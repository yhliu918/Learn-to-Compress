
#ifndef PIECEWISE_H_
#define PIECEWISE_H_

#include "common.h"
#include "codecs.h"
#include "time.h"
#include "bit_read.h"
#include "bit_write.h"
#include "caltime.h"
#define INF 0x7f7fffff

namespace Codecset {






double cal_score_tmp(uint32_t* a, int N){
  int * b = (int*)malloc(N * sizeof(int));
  uint32_t amax = a[N-1];
  uint32_t amin = a[N-1];
  for(int i=0;i<N-1;i++){
    if(amax<a[i]){
      amax = a[i];
    }
    if(amin>a[i]){
      amin = a[i];
    }
    b[i] = a[i+1]-a[i];
  }
  int best_delta = (amax-amin)/(N-1);
  long long sum =0;
  for(int i=0;i<N-1;i++){
    sum += abs(abs(b[i]) - best_delta);
  }
  if(amax -amin == 0){
    return 1;
  }
  double score = (double)sum / ((amax-amin));
  score/=(double)1;
  free(b);
  //std::cout<<"score "<<score<<std::endl;
  return 1-score;
}
class piecewise : public IntegerCODEC {
public:
  using IntegerCODEC::encodeArray;
  using IntegerCODEC::decodeArray;
  using IntegerCODEC::randomdecodeArray;
  using IntegerCODEC::encodeArray8;
  using IntegerCODEC::decodeArray8;
  using IntegerCODEC::randomdecodeArray8;
  using IntegerCODEC::init;

  std::vector<uint8_t*> block_start_vec;
  std::vector<uint32_t> segment_index;
  uint32_t counter = 0;
  uint32_t total_byte =0;
  int maxerror =(1U<<10)-1;
  int block_num;
  int block_size;

//start_index + bit + theta0 + theta1 + numbers + delta
void init( int blocks,  int blocksize,int delta){
      block_num=blocks;
      block_size=blocksize;
      maxerror = delta;

}

float epsilon = 0.;
bool gt(float a, float b){
// greater than
    return (a > b + epsilon);
}
uint32_t lower_bound( uint32_t v,uint32_t len)
{
    uint32_t m;
    uint32_t x=0;
    uint32_t y=len-1;
    while(x <= y)
    {

        m = x+(y-x)/2;
        if(v<segment_index[m]) y = m-1;
        else x = m+1;
    }
    return y;

}

uint8_t * encodeArray8(uint32_t *in, const size_t length,uint8_t *res, size_t nvalue) {

    std::vector<uint32_t> indexes;
    for(uint32_t i = 0; i < nvalue; i++){
        indexes.push_back(i);
    }
    double high_slope = (double)INF;
    double low_slope = 0.;
    long long origin_key = in[0];
    int origin_index = indexes[0];
    int end_index = indexes[0];
    int total_index =0;
    for (int i = 1; i < (long long)nvalue; i++){
        long long key = in[i];
        int id = indexes[i];
        double tmp_point_slope = ((key - origin_key)+0.0) / ((id - origin_index)+0.0);



        if (tmp_point_slope >= low_slope && tmp_point_slope <= high_slope  ){
        //if (gt(tmp_point_slope, low_slope) && gt(high_slope, tmp_point_slope)){
            double tmp_high_slope = ((key + maxerror - origin_key)+0.0) / ((id - origin_index)+0.0);
            double tmp_low_slope = ((key - maxerror - origin_key)+0.0) /((id - origin_index)+0.0);
            if (tmp_low_slope<0.0){
                //std::cout<<"zero!"<<std::endl;
                tmp_low_slope=0.0;
            }
            //std::cout<<tmp_high_slope<<","<<tmp_low_slope<<std::endl;

           if(tmp_high_slope<=high_slope){
           //if(gt(high_slope, tmp_high_slope)){
                high_slope = tmp_high_slope;
            }
            if(low_slope<=tmp_low_slope){
            //if(gt(tmp_low_slope, low_slope)){
                low_slope = tmp_low_slope;
            }
            end_index = id;

        }
        else{

            double slope = (high_slope + low_slope) / 2.;
            int max_error = 0;

            if (end_index == origin_index){
                slope = 1.;
            }
            int seg_len = end_index - origin_index + 1;
            int * delta = static_cast<int *>(malloc(seg_len * sizeof(int)));


            for (int j = origin_index; j <= end_index; j++ ){
                long long  predict = (long long)(in[origin_index]+slope*(double)(j-origin_index));
                long long truth = (long long)in[j];
                int tmp = abs(predict-truth);

                delta[j-origin_index] = (int) truth-predict;
                total_index++;

                if (tmp > max_error){
                    max_error = tmp;

                }
            }
            
            //std::cout<<cal_score_tmp(in+origin_index,(end_index - origin_index + 1))<<std::endl;
            int tmp_bit = 0;
            if(max_error){
                tmp_bit = bits(max_error)+1;
            }

            // std::cout<<tmp_bit<<" ("<<origin_index<<" , "<<end_index<<") "<<(end_index - origin_index + 1)<<" "<<in[origin_index]<<" "<<in[end_index]<<std::endl;
            uint8_t * descriptor = (uint8_t*)malloc((end_index - origin_index + 1) * sizeof(uint64_t)+30);

            uint8_t *out = descriptor;

            memcpy(out,&origin_index,sizeof(origin_index));
            out+=sizeof(origin_index);

            out[0]=(uint8_t)tmp_bit;
            out++;

            uint32_t theta0 =  (long long)in[origin_index] ;
            uint32_t numbers = (end_index - origin_index + 1);
            memcpy(out,&theta0,sizeof(theta0));

            out+=sizeof(theta0);
            memcpy(out,&slope,sizeof(slope));
            out+=sizeof(slope);
            

            //std::cout<<"bit_length: "<<tmp_bit<<" start: "<<origin_index<<" end: "<<end_index<<" slope: "<<slope<<std::endl;
            if(tmp_bit){
                out = write_delta_T(delta, out, tmp_bit, seg_len);
            }
            
            free(delta);

            if(numbers<=2){
                counter++;
                // std::cout<<out-descriptor<<std::endl;
            }
            descriptor=(uint8_t*)realloc(descriptor, (out-descriptor));
            block_start_vec.push_back(descriptor);
            segment_index.push_back(origin_index);

            total_byte +=(out - descriptor);
            //std::cout<<"total_byte"<<total_byte<<std::endl;

            //std::cout<<"tmp_bit"<<unsigned(tmp_bit)<<" numbers "<<numbers<<" slope "<<slope<<std::endl;
            //std::cout<<"bit_length: "<<tmp_bit<<" start: "<<origin_index<<" end: "<<end_index<<" slope: "<<slope<<std::endl;
            high_slope = (double)INF;
            low_slope = 0.0;
            origin_index = id;
            origin_key = key;
            end_index = id;


        }

    }

    double slope = (high_slope + low_slope) / 2.;
    if (end_index == origin_index){
        slope = 1.;
    }

    int seg_len = end_index - origin_index + 1;
    int * delta = static_cast<int *>(malloc(seg_len * sizeof(int)));
    int max_error = 0;
    for (int j = origin_index; j <= end_index; j++ ){
        long long  predict = (long long)in[origin_index] + (long long)(slope*(double)(j-origin_index));
        long long truth = (long long)in[j];
        int tmp = abs(predict-truth);
        delta[j-origin_index] = (int) truth-predict;

        total_index++;
        if (tmp > max_error){
            max_error = tmp;

        }
    }

    uint8_t tmp_bit = 0;
    if(max_error){
        tmp_bit = bits(max_error)+1;
    }


    uint8_t * descriptor = (uint8_t*)malloc((end_index - origin_index + 1) * sizeof(uint64_t)+30);
    uint8_t *out = descriptor;
    uint32_t start_ind = origin_index;
    memcpy(out,&start_ind,sizeof(start_ind));
    out+=sizeof(start_ind);

    out[0]=(uint8_t)tmp_bit;
    out++;

    uint32_t theta0 = (long long)in[origin_index] ;
    uint32_t numbers = (end_index - origin_index + 1);
    memcpy(out,&theta0,sizeof(theta0));
    out+=sizeof(theta0);
    memcpy(out,&slope,sizeof(slope));
    out+=sizeof(slope);

    if(tmp_bit){
        out = write_delta_T(delta, out, tmp_bit, seg_len);
    }
    free(delta);

    descriptor=(uint8_t*)realloc(descriptor, (out-descriptor));
    block_start_vec.push_back(descriptor);
    segment_index.push_back(start_ind);
    total_byte +=(out - descriptor);

    // std::cout<<"small frag number "<< counter<<std::endl;

    return res;

}

uint32_t *decodeArray8( uint8_t *in, const size_t length, uint32_t *out, size_t nvalue) {
//start_index + bit + theta0 + theta1 + numbers + delta
    uint32_t * tmpout = out;
    for(int i=0;i<(int)block_start_vec.size();i++){
        uint8_t * this_block = block_start_vec[i];
        uint8_t * tmpin=this_block;
        uint32_t theta0;
        double theta1;
        uint8_t maxerror;
        uint32_t start_ind;
        uint32_t numbers;

        memcpy(&start_ind,tmpin,4);
        tmpin+=4;
        memcpy(&maxerror,tmpin,1);
        tmpin++;
        memcpy(&theta0,tmpin,4);
        tmpin+=4;
        memcpy(&theta1,tmpin,sizeof(double));
        tmpin+=sizeof(double);
        if(i==block_start_vec.size()-1){
            numbers = length - segment_index[i];
        }
        else{
            numbers = segment_index[i+1] - segment_index[i];
        }

        /*
        if(start_ind ==232){
        std::cout<<"start_ind "<<start_ind<<" maxerror "<<unsigned(maxerror)<<" theta0 "<<theta0<<" theta1 "<<theta1<<" numbers "<<numbers<<std::endl;
        }
        */
        if(numbers ==1){
            tmpout[0]=theta0;
            tmpout++;


        }
        else{
            tmpout=read_all_bit(tmpin ,0,start_ind,numbers,maxerror,theta1,theta0, tmpout);
        }
        
    }





    return out;
}

uint32_t randomdecodeArray8( uint8_t *in, const size_t l, uint32_t *out, size_t nvalue){

    uint32_t length=segment_index.size();
    uint8_t * this_block = block_start_vec[lower_bound(l,length)];

    uint8_t * tmpin = this_block;
    uint32_t theta0;
    double theta1;
    uint8_t maxerror;
    uint32_t start_ind;
    uint32_t numbers;

    memcpy(&start_ind,tmpin,4);
    tmpin+=4;
    memcpy(&maxerror,tmpin,1);
    tmpin++;
    memcpy(&theta0,tmpin,4);
    tmpin+=4;
    memcpy(&theta1,tmpin,sizeof(theta1));
    tmpin+=sizeof(theta1);

    //std::cout<< "indexing "<<l<<std::endl;
    uint32_t tmp = 0;
    if(maxerror){
        tmp = read_bit_fix_T(tmpin ,maxerror, l-start_ind, theta1,theta0, 0);
    }
    else{
        tmp = (long long)(theta0+theta1*(double)(l-start_ind));
    }
    
    

    return tmp;

}
uint64_t summation( uint8_t *in, const size_t l, size_t nvalue){

    return 0;
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
    std::cout<<"Total block num is "<<block_start_vec.size()<<std::endl;
     return total_byte;
 }

void destroy(){
for(int i=0;i<(int)block_start_vec.size();i++){
       free(block_start_vec[i]);
   }

}
std::string name() const {
    return "piecewise";
}

};



} // namespace FastPFor

#endif /* SIMDFASTPFOR_H_ */
