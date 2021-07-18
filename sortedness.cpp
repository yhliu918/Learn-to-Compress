// sortedness v.s. compression rate

#include "../headers/common.h"
#include "../headers/codecfactory.h"
#include "../headers/caltime.h"
#include "../headers/lr.h"

using namespace Codecset;

int random(int m){
        return rand()%m;
}

long long count =0;
void merge_solve(uint32_t* A, int start, int end){
    if(start >= end) return;                                        //只剩下一个元素时

    int s = start;
    int e = end;
    int mid = (start + end) / 2;                                    //中间的元素

    merge_solve(A, s, mid);
    merge_solve(A, mid + 1, e);

    //回溯阶段,两个有序数组合并，需要一个新的数组记住拍好的序之后的值
    uint32_t* tmp = new uint32_t[end - start + 1];                            //存临时结果的数组
    int j = 0;
    int s2 = mid + 1;

    while(s <= mid && s2 <= e){
        if(A[s] <= A[s2])
            tmp[j++] = A[s++];
        else{
            count += (mid - s + 1);                                 //逆序对数等于前半中剩余元素的个数
            tmp[j++] = A[s2++];
        }
    }

    while(s <= mid)
        tmp[j++] = A[s++];
    while(s2 <= e)
        tmp[j++] = A[s2++];

    //复制数据
    s = start;
    while(s <= e){
        A[s] = tmp[s - start];
        s++;
    }
    delete[] tmp;
}

double cal_inverse(uint32_t a[], int N){
    count =0;
    uint32_t* tmp = new uint32_t[N];
    for(int i=0;i<N;i++){
      tmp[i]=a[i];
    }
    merge_solve(tmp, 0, N);
    
    double score =(double)count/(double)N;
    score =score/(double)(N-1);
    score*=2;
    delete [] tmp;
    return (fabs(1-2.0*score)*fabs(1-2.0*score));   // inverse = C(n,2) ==1 inverse == 0 -->0 0 random --> 0.5
}

double cal_score(uint32_t a[], int N){
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
    return 0;
  }
  double score = (double)sum / ((amax-amin));
  score/=(double)(2*pow(10,5));
  free(b);
  //std::cout<<"score "<<score<<std::endl;
  return 1-score;
}

IntegerCODEC &codec1 = *CODECFactory::getFromName("piecewise_fix");
IntegerCODEC &codec2 = *CODECFactory::getFromName("FOR");

uint32_t * data =NULL;
int N = 2000000;
  int blocks =1000;
  int block_size = N/blocks;
int partition(uint32_t a[], int l, int r)
{
    int i = l;
    int j = r + 1;
    int v = a[l];
    
    while (true) {
        while (a[--j] > v);
        while (i < j && a[++i] < v);
        
        if (i == j)
        {
            break;
        }
        std::swap(a[i], a[j]);
    }
    
    std::swap(a[i], a[l]);
    
    return i;
}
void compress (uint32_t a[],int N){
  int totalsize = 0;
  for(int i=0;i<blocks;i++){
    uint8_t * descriptor = (uint8_t*)malloc(block_size * sizeof(uint64_t));
    uint8_t * res = descriptor;
    res = codec1.encodeArray8(data+(i*block_size),block_size ,descriptor,i);
    totalsize += (res-descriptor);
    free(descriptor);
  }
  double compressrate = (totalsize)*100.0  / (4*N*1.0);
  std::cout << std::setprecision(8)<< compressrate << std::endl;

  totalsize = 0;
  for(int i=0;i<blocks;i++){
    uint8_t * descriptor = (uint8_t*)malloc(block_size * sizeof(uint64_t));
    uint8_t * res = descriptor;
    res = codec2.encodeArray8(data+(i*block_size),block_size ,descriptor,i);
    totalsize += (res-descriptor);
    free(descriptor);
  }
  compressrate = (totalsize)*100.0  / (4*N*1.0);
  std::cout << std::setprecision(8)<< compressrate << std::endl;
  
}
int globalk =0;
double last_arl =0;
void quick_sort(uint32_t a[], int l, int r,int N,double extent)
{
  
  //if(random(10) || last_arl>0.9){
    double arl =cal_inverse(a,N);
    if(fabs(arl - last_arl)>0.0001){
      std::cout<<arl<<std::endl;
      last_arl = arl;
      compress(a,N);
    }
  //}
  
  /*
  if(random(10000) == 0){
    double arl = cal_score(a,N);
    std::cout<<arl<<std::endl;

  }
  */
  
  globalk++;
  //if(arl>extent){
    //  return;
  //}
  //compress(a,N);
    if (l < r) {
      
      
        int m = partition(a, l, r);
        
        quick_sort(a, l, m - 1,N,extent);
        quick_sort(a, m + 1, r,N,extent);
    }
}


int main() {
  
  // We pick a CODEC
  
  data = (uint32_t*)malloc(N * sizeof(uint32_t));
  //std::ifstream srcFile("../data/sortedness_normal/blocksort_1000000_data_200M.txt",std::ios::in); 
  std::ifstream srcFile("/root/database/uwg/arrow/data/build/data_200M.txt",std::ios::in); 
  if(!srcFile) { 
      std::cout << "error opening source file." << std::endl;
      return 0;
  }
  for(int i=0;i<N;i++){
    uint32_t num;
    srcFile>>num;
    data[i]=num;
  }
  srcFile.close();
  

  
  codec1.init(blocks,block_size,32);
  codec2.init(blocks,block_size,32);

  quick_sort(data, 0, N - 1,N,0);
  //std::cout<<globalk<<std::endl;
  



  
}
