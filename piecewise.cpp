#include <stdio.h>
#include<iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <cmath>
#include <time.h>
#include<algorithm>
#include<queue>
#include<vector>
#include <getopt.h>
#include <assert.h>
#include <time.h>
#include <sys/time.h>
#include "bit_opt.h"
#define INF 0x7f7fffff

uint32_t NUM_KEYS=200000000;
std::vector<uint32_t> segindex; 
double getNow() {
  struct timeval tv;
  gettimeofday(&tv, 0);
  return tv.tv_sec + tv.tv_usec / 1000000.0;
}
// we can calculate the first_byte via start index (each segment's delta would take a interger times of bytes byte)
struct Segment{
    uint32_t start_index;
    float slope;
    uint8_t max_error; // how many bits we need to save this segment
    uint32_t start_byte;
};

uint32_t lower_bound( uint32_t v)
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

class Piecewise
{
    std::vector<uint32_t> keys;
    std::vector<uint32_t> indexes;
    std::vector<Segment> q;
    int maxerror ;
    uint32_t totalbit=0;
    uint8_t * compressed = static_cast<uint8_t *>(malloc(NUM_KEYS * sizeof(uint32_t)));
    int * delta_origin = static_cast<int *>(malloc(NUM_KEYS * sizeof(int)));
    int * delta_after = static_cast<int *>(malloc(NUM_KEYS * sizeof(int)));
    //uint8_t * compressed = malloc(NUM_KEYS * sizeof(uint32_t));
    //int * delta_origin = malloc(NUM_KEYS * sizeof(int));
    
    
    public:
    
    uint32_t getTotalbit();
    void setMaxError(int error);
    void setkeys(char* s);
    void setindex();
    void fit();
    void printSeg();
    void decompress();
    void ra_decompress();
    void decompress_all();
    int write_delta(int *in, int start_byte, uint8_t l, int numbers);
    int seglen();
    void printCompressed();
};
void Piecewise::printCompressed(){
    for(int i =0; i <NUM_KEYS;i++){
        std::cout<<unsigned(compressed[i])<<std::endl;
    }
}
uint32_t Piecewise::getTotalbit(){
    return totalbit;
}

void Piecewise::setMaxError(int error){
    maxerror = error;
    std::cout<<maxerror<<std::endl;
}

void Piecewise::setkeys(char* s){
    char filename[50]={0};
    char head[10]="data/";
    char tail[20]="_200M_uint32.txt";
    int k = 0;
    for(int i =0; i< strlen(head);i++){
        filename[k]=head[i];
        k++;
    }
    for(int i =0; i< strlen(s);i++){
        filename[k]=s[i];
        k++;
    }
    for(int i =0; i< strlen(tail);i++){
        filename[k]=tail[i];
        k++;
    }

    std::ifstream infile(filename);
    for(int i = 0; i < NUM_KEYS; i++){
        uint32_t tmpnum;
        infile>>tmpnum;
        keys.push_back(tmpnum);
        //std::cout<<keys[i]<<std::endl;
        
    } 
}

void Piecewise::setindex(){
    for(uint32_t i = 0; i < NUM_KEYS; i++){
        indexes.push_back(i);
    }   
}


void Piecewise::printSeg(){
    for(int i=0;i<q.size();i++){
        Segment tmpseg = q[i];
        printf(" %lu, %lu, %.1f, %lu\n ",tmpseg.start_index, tmpseg.start_byte,  tmpseg.slope, tmpseg.max_error);

    }

}
int Piecewise::seglen(){
    return q.size();
}

float epsilon = 0.;
bool gt(float a, float b){
// greater than
    return (a > b + epsilon);
}

int Piecewise::write_delta(int *in,int start_byte, uint8_t l, int numbers){
    int code =0;
    int occupy = 0;
    int write_byte=start_byte;
    uint8_t left_val = 0;

    for(int i = 0; i < numbers;i++){

        code = left_val;
        while(occupy<8){
            
            bool sign = 1;
            if (in[i] <= 0){
                sign = 0;
                in[i]=-in[i];
            }
            //std::cout<<(in[i] & ((1U<<(l-1))-1))<<","<<(sign<<(l-1))<<","<<((in[i] & ((1U<<(l-1))-1)) + (sign<<(l-1)))<<std::endl;
            /*
            if(start_byte==143750518){
                std::cout<<"befor code:"<<unsigned((uint8_t)code)<<std::endl;
            }
            */
            uint8_t value1= ((in[i] & ((1U<<(l-1))-1)) + (sign<<(l-1)));
            code += (value1<<occupy);
            occupy += l;
            /*
            if(start_byte==143750518){
                std::cout<<i<<","<<occupy<<", code:"<<unsigned((uint8_t)code)<<","<<in[i]<<","<<unsigned(value1)<<","<<sign<<std::endl;
            }
            */
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
        
        compressed[write_byte]= unsigned((uint8_t)code);
        write_byte++;
        if(i==numbers){
            if (left_val){
                compressed[write_byte]=(uint8_t)left_val;
                write_byte++;
            }
            if(left_val==0&&occupy!=0){
                compressed[write_byte]=(uint8_t)left_val;
                write_byte++;
            }
            
        }
        i--;
    }
    /*
    if(start_byte==143750518){
        for(int i=start_byte;i<write_byte;i++){
            std::cout<<"write_byte: "<<i<<","<<unsigned(compressed[i])<<std::endl;
        }
    }
    */
    
    return write_byte;

}
void Piecewise::fit(){
    float high_slope = (float)INF;
    float low_slope = 0.;
    long long origin_key = keys[0];
    int origin_index = indexes[0];
    int end_index = indexes[0];
    int start_byte=0;
    int total_index =0;
    for (int i = 1; i < NUM_KEYS; i++){
        long long key = keys[i];
        int id = indexes[i];
        float tmp_point_slope = ((key - origin_key)+0.0) / ((id - origin_index)+0.0);
        
//        if (tmp_point_slope >= low_slope && tmp_point_slope <= high_slope){
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
            uint8_t max_error = 0;

            if (end_index == origin_index){
                slope = 1.;
            }
            int seg_len = end_index - origin_index + 1;
            int * delta = static_cast<int *>(malloc(seg_len * sizeof(int)));
            
            max_error = 0;
            for (int j = origin_index; j <= end_index; j++ ){
                long long  predict = (long long)keys[origin_index] + (long long)(slope*(float)(j-origin_index));
                long long truth = (long long)keys[j];
                int tmp = abs(predict-truth);
                
                delta[j-origin_index] = (int) predict-truth;
                delta_origin[total_index]=delta[j-origin_index];
                total_index++;
                // std::cout<<"index:"<<j<<"key:"<<(long long)truth<<std::endl;
                /*
                if (tmp > maxerror){
                    std::cout<<"predict is "<<predict<<std::endl;
                    std::cout<<"truth is "<<truth<<std::endl;
                    std::cout<<"tmp is "<<tmp<<std::endl;
                    std::cout<<"slope is "<<(double)slope<<std::endl;
                    std::cout<<"startkey is "<<keys[origin_index]<<std::endl;
                    std::cout<<"start index is "<<origin_index<<std::endl;
                    std::cout<<"end index is "<<end_index<<std::endl;
                    std::cout<<"this index is "<<j<<std::endl;
                    std::cout<<"Something wrong!"<<std::endl;
                    
                    std::cout<<std::endl<<"tmp point slope:"<<tmp_point_slope<<std::endl;
                    std::cout<<"low slope:"<<(double)low_slope<<std::endl;
                    std::cout<<"high slope:"<<(double)high_slope<<std::endl;
                    std::cout<<"diff:"<<high_slope-low_slope<<std::endl;
                    for(int iii=origin_index;iii<=end_index;iii++){
                        std::cout<<"index:"<<iii<<"key:"<<(long long)keys[iii]<<std::endl;
                    } 
                    
                    std::cout<<(keys[j]+0.0)<<" versus "<<long(keys[j])<<std::endl;
                    exit(0);
                }
                */
                if (tmp > max_error){
                    max_error = tmp;
                    
                }
            }
            if (max_error > 0.01){
                totalbit += (ceil(log2(max_error))+2)*(end_index-origin_index+1);//have a bit of sign
            }
            else{
                totalbit += 2 * (end_index - origin_index + 1);
            }
            totalbit +=104 ;
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
            /*
            if(start_byte== 199999998){
                for(int j = 0; j < (end_index - origin_index + 1); j++ ){
                    std::cout<<delta[j]<<std::endl;
                }
            }
            */
            //std::cout<<start_byte<<","<<(end_index - origin_index + 1)<<std::endl;
            start_byte = write_delta(delta,start_byte, tmp_bit, (end_index - origin_index + 1));
            

            
            q.push_back(newseg);
            segindex.push_back(origin_index);
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
    uint8_t max_error = 0;
    for (int j = origin_index; j <= end_index; j++ ){
        long long  predict = (long long)keys[origin_index] + (long long)(slope*(float)(j-origin_index));
        long long truth = (long long)keys[j];
        int tmp = abs(predict-truth);     
        delta[j-origin_index] = (int) predict-truth;
        delta_origin[total_index]=delta[j-origin_index];
        total_index++;
        if (tmp > max_error){
            max_error = tmp;   
        }
    }
    if (max_error > 0.01){
        totalbit += (ceil(log2(max_error))+2)*(end_index-origin_index+1);//have a bit of sign
    }
    else{
        totalbit += 2 * (end_index - origin_index + 1);
    }
    totalbit +=104 ;
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
    

    //std::cout<<start_byte<<","<<(end_index - origin_index + 1)<<std::endl;
    start_byte = write_delta(delta,start_byte, tmp_bit, (end_index - origin_index + 1));
    q.push_back(newseg);
    segindex.push_back(origin_index);
    free(delta);
    
}

void Piecewise::ra_decompress(){
    std::cout<<"random access decompressing!"<<std::endl;
    double totaltime=0;
    for(int i=0;i<NUM_KEYS;i++){
      
      double start_time = getNow();
      
      Segment tmpseg = q[lower_bound(i)];

      delta_after[i]=read_bit(compressed ,tmpseg.max_error , i-tmpseg.start_index, tmpseg.start_byte);
      double end_time = getNow();
      totaltime += (end_time - start_time);
      if(delta_after[i]!=delta_origin[i]){
        std::cout<<lower_bound(i)<<std::endl;
        std::cout<<tmpseg.start_byte<<"start is: "<<tmpseg.start_index<<","<<i<<","<<delta_after[i]<<","<<delta_origin[i]<<std::endl;
        std::cout<<"something wrong! decompress failed!"<<std::endl;
        return;
      }
    }

    //assert(delta_origin == delta_after);
    std::cout<<"total decompress time:"<<totaltime<<std::endl;
    std::cout<<"average decompress time:"<<(double)(totaltime/NUM_KEYS)*1000000000<<std::endl;
    std::cout<<" decompress speed:"<<(double)(NUM_KEYS/(totaltime*1000))<<std::endl;
    free(delta_origin);
    free(delta_after);

}
void Piecewise::decompress_all(){
    std::cout<<"decompressing all!"<<std::endl;
    int start_byte=0;
    int end_index =0;
    int start_index = 0;
    double totaltime=0;
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
            end_index = NUM_KEYS-1;
        }
        double start_time = getNow();
        read_all_bit(compressed ,start_byte,start_index, (end_index-start_index+1),tmpseg.max_error, delta_after);
        double end_time = getNow();
        totaltime += (end_time - start_time);
        for (int i=start_index;i<end_index;i++){
            
            if(delta_after[i]!=delta_origin[i]){
                std::cout<<start_byte<<"start is: "<<start_index<<", end is: "<<end_index<<","<<i<<","<<delta_after[i]<<","<<delta_origin[i]<<std::endl;
                std::cout<<"something wrong! decompress failed!"<<std::endl;
                return;
            }
        }
        
    }

    //assert(delta_origin == delta_after);
    std::cout<<"total decompress time:"<<totaltime<<std::endl;
    std::cout<<"average decompress time:"<<(double)(totaltime/NUM_KEYS)*1000000000<<std::endl;
    std::cout<<" decompress speed:"<<(double)(NUM_KEYS/(totaltime*1000))<<std::endl;
    free(delta_origin);
    free(delta_after);

}
void Piecewise::decompress(){
    std::cout<<"decompressing!"<<std::endl;
    int start_byte=0;
    int start_index = 0;
    int end_index =0;
    int total_index=0;
    double totaltime=0;
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
            end_index = NUM_KEYS-1;
        }
        for (int i=0;i<end_index-start_index+1;i++){
            double start_time = getNow();
            delta_after[total_index]=read_bit(compressed ,tmpseg.max_error , i , start_byte);
            double end_time = getNow();
            totaltime += (end_time - start_time);
            if(delta_after[total_index]!=delta_origin[total_index]){
                std::cout<<start_byte<<"start is: "<<start_index<<", end is: "<<end_index<<","<<total_index<<","<<delta_after[total_index]<<","<<delta_origin[total_index]<<std::endl;
                std::cout<<"something wrong! decompress failed!"<<std::endl;
                return;
            }
            total_index++;
        }
    }

    //assert(delta_origin == delta_after);
    std::cout<<"total decompress time:"<<totaltime<<std::endl;
    std::cout<<"average decompress time:"<<(double)(totaltime/NUM_KEYS)*1000000000<<std::endl;
    std::cout<<" decompress speed:"<<(double)(NUM_KEYS/(totaltime*1000))<<std::endl;
    free(delta_origin);
    free(delta_after);

}
int char2int(const char* str) {
    const char* p = str;
    bool neg = false;
    int res = 0;
    if (*str == '-' || *str == '+') {
        str++;
    }

    while (*str != 0) {
        if (*str < '0' || *str > '9') {
            break;
        }
        res = res * 10 + *str - '0';
        str++;
    }

    if (*p == '-') {
        res = -res;
    }
    return res;
}



int main(int argc, char *argv[])
{
    
    static struct option long_options[] =
    {
        {"file", required_argument,NULL,'f'},
        {"error",required_argument,NULL,'e'},
    };
    
    Piecewise model;
    int opt_index = 0;
    int c; 
    bool defaultkeyflag=true;
    bool defaulterrorflag=true;
    while(1)
    {
        int opt_index = 0;
        c = getopt_long(argc, argv,"f:e:", long_options,&opt_index);

        if(-1 == c)
        {
            break;
        }
        
        switch(c)
        {
            case 'f':
                std::cout<<"Preparing data..."<<std::endl;
                model.setkeys(optarg);
                defaultkeyflag=false;
                break;
            case 'e':
                model.setMaxError(char2int(optarg));
                defaulterrorflag=false;
                break;
            default:
                std::cout << "???" << std::endl;
                break;
        }
    }
    if(defaultkeyflag){
        std::cout<<"Preparing data..."<<std::endl;
        model.setkeys("normal");
    }
    if(defaulterrorflag){
        model.setMaxError(127);
    }
    
    
    model.setindex();
    std::cout<<"Model Training..."<<std::endl;
    model.fit();
    int blocks = model.seglen();
    //model.printSeg();
    uint32_t total = model.getTotalbit();
    // model.decompress();
    //model.decompress_all();
    model.ra_decompress();
    //model.printCompressed();
    
    std::cout<<"Splitting into "<<blocks<<" pieces"<<std::endl;
    std::cout<<"compressing "<<NUM_KEYS<<"itegers ("<<NUM_KEYS*4<<" bytes) into "<<total/8<<" bytes, compression rate is "<<(float)total/(8.0*NUM_KEYS*4)<<std::endl;

    return 0;
}
