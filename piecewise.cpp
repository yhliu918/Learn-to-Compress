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

#define INF 0x7f7fffff

uint32_t NUM_KEYS=200000000;

struct Segment{
    uint32_t start;
    uint32_t end;
    double slope;
    int max_error;

};

class Piecewise
{
    std::vector<uint32_t> keys;
    std::vector<uint32_t> indexes;
    std::queue<Segment> q;
    int maxerror ;
    uint32_t totalbit=0;
    
    public:
    
    uint32_t getTotalbit();
    void setMaxError(int error);
    void setkeys(char* s);
    void setindex();
    void fit();
    void printSeg();
    int seglen();
};

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

    while(!q.empty()){
        Segment tmpseg = q.front();
        printf(" %lu, %lu, %.1f, %lu\n ",tmpseg.start, tmpseg.end, tmpseg.slope, tmpseg.max_error);
        q.pop();
    }
}
int Piecewise::seglen(){
    return q.size();
}
void Piecewise::fit(){
    double high_slope = (double)INF;
    double low_slope = 0.;
    uint32_t origin_key = keys[0];
    uint32_t origin_index = indexes[0];
    uint32_t end_index = indexes[0];
    for (uint32_t i = 1; i < NUM_KEYS; i++){
        uint32_t key = keys[i];
        uint32_t id = indexes[i];
        double tmp_point_slope = (float)(int)(key - origin_key) / (float)(int)(id - origin_index);
        
        if (tmp_point_slope >= low_slope && tmp_point_slope <= high_slope){
            double tmp_high_slope = (double)(int)((key + maxerror) - origin_key) / (double)(int)(id - origin_index);
            double tmp_low_slope = (double)(int)((key - maxerror) - origin_key) /(double) (int)(id - origin_index);
            
            //std::cout<<tmp_low_slope<<","<<tmp_high_slope<<std::endl;
            if(tmp_high_slope<high_slope){
                high_slope = tmp_high_slope;
            }
            if(low_slope<tmp_low_slope){
                low_slope = tmp_low_slope;
            }
            
            end_index = id;
            
        }
        else{
            
            double slope = (high_slope + low_slope) / 2.;
            
            if (end_index == origin_index){
                slope = 1;
            }
            int max_error = 0;
            for (uint32_t j = origin_index+1; j <= end_index; j++ ){
                int tmp = abs((double)keys[origin_index] + (slope*(j-origin_index))-(double)keys[j]);
                //std::cout<<"predict is "<<keys[origin_index] + (int)(slope*(j-origin_index))<<std::endl;
                //std::cout<<"truth is "<<keys[j]<<std::endl;
                //if (tmp > maxerror){
                    //std::cout<<"Something wrong!"<<std::endl;
                //}
                if (tmp > max_error){
                    max_error = tmp;
                    
                }
            }
            if (max_error > 0.01){
                totalbit += (ceil(log2(max_error))+1)*(end_index-origin_index+1);
            }
            else{
                totalbit += (end_index - origin_index + 1);
            }
            totalbit += 160;
            Segment newseg;
            newseg.start = origin_index;
            newseg.end = end_index;
            newseg.slope = slope;
            newseg.max_error = max_error;
            q.push(newseg);
            
            high_slope = (double)INF;
            low_slope = 0.;
            origin_index = id;
            origin_key = key;
            end_index = id;
        }

    }
    
    double slope = (high_slope + low_slope) / 2.;
    if (end_index == origin_index){
        slope = 1;
    }

    int max_error = 0;
    for (uint32_t j = origin_index+1; j <= end_index; j++ ){
        int tmp = abs((double)keys[origin_index] + (slope*(j-origin_index))-(double)keys[j]);
        if (tmp > max_error){
        
            max_error = tmp;
            //std::cout<<j<<","<<max_error<<std::endl;
        }
    }
    if (max_error > 0.01){
        totalbit += (ceil(log2(max_error))+1)*(end_index-origin_index+1);
    }
    else{
        totalbit += (end_index - origin_index + 1);
    }
    totalbit += 160;
    Segment newseg;
    newseg.start = origin_index;
    newseg.end = end_index;
    newseg.slope = slope;
    newseg.max_error = max_error;
    q.push(newseg);
    
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
    std::cout<<"Splitting into "<<blocks<<" pieces"<<std::endl;
    std::cout<<"compressing "<<NUM_KEYS<<"itegers ("<<NUM_KEYS*4<<" bytes) into "<<total/8<<" bytes, compression rate is "<<(float)total/(8.0*NUM_KEYS*4)<<std::endl;

    return 0;
}
