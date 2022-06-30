#include "../headers/RANSAC.h"
#include <cstring>
#include <math.h>
#include <algorithm>

RANSAC::RANSAC()
  :probability_(0.99), max_iterations_ (10)
{

}

int RANSAC::compute(lr& mylr, double x[], double y[],int numbers,int num_for_estimate,bool *choose)
{

  int num_data = numbers;

  int iterations = 0;  
  double k = 100;  


  double log_probability  = log (1.0 - probability_);
  double one_over_indices = 1.0 / static_cast<double> (num_data);
  
  

  short *best_votes = new short[num_data];
  short *cur_votes = new short[num_data]; 

  SubSetIndexComparator sub_set_index_comparator(num_for_estimate);
  std::set<int *, SubSetIndexComparator > chosen_sub_sets(sub_set_index_comparator);  


  int* cur_inti_sub_set_indexs = NULL;
  

  int best_model_num = -1;  
  int maybe_inliers_num = 0;  
  std::vector<double> maybe_model;  

  std::vector<int> shuffled_indices(num_data);

  while (iterations < k)
  {
      
      maybe_inliers_num = 0;
      cur_inti_sub_set_indexs = new int[num_for_estimate];
          

      maybe_model.clear();
    
      for (int i=0;i < (int)shuffled_indices.size(); i++)
      {
        shuffled_indices[i]=i;
      }

      int max_index = num_data-1;
      for (int i=0; i<num_for_estimate; i++)
      {
        std::swap(shuffled_indices[i], shuffled_indices[i + rand() % (numbers - i)]);
      }
    
      memset(cur_votes, 0, num_data*sizeof(short));

      for (int i=0; i<num_for_estimate; i++)
      {
        cur_inti_sub_set_indexs[i] = shuffled_indices[i];
        cur_votes[shuffled_indices[i]] = 1;
      }
      maybe_inliers_num = num_for_estimate;


      std::pair< std::set<int *, SubSetIndexComparator >::iterator, bool > res = chosen_sub_sets.insert(cur_inti_sub_set_indexs);

      if (res.second)
      {

        std::vector<double> tmpx;
        std::vector<double> tmpy;
        for (int i=0; i<num_for_estimate; i++)
        {
          tmpx.emplace_back(x[cur_inti_sub_set_indexs[i]]);
          tmpy.emplace_back(y[cur_inti_sub_set_indexs[i]]);
        }


        mylr.caltheta(tmpx,tmpy,num_for_estimate);

        for(int i=0; i<num_data; i++)
        {
          if(0 == cur_votes[i] && 
            mylr.agree(x[i],y[i]))
          {
            cur_votes[i] = 1;
            maybe_inliers_num++;          
          }
        }

        if (maybe_inliers_num > best_model_num)
        {
          best_model_num = maybe_inliers_num;
          memcpy(best_votes, cur_votes, num_data*sizeof(short));
        }

      //k=log(1-p)/log(1-pow(w,n))
        double w = static_cast<double> (best_model_num) * one_over_indices;
        double p_no_outliers = 1.0 - std::pow(w, static_cast<double> (maybe_inliers_num));
        if(abs(p_no_outliers-1.0)<std::numeric_limits<double>::epsilon ()){
          break;
        }
        p_no_outliers = (std::max) (std::numeric_limits<double>::epsilon (), p_no_outliers);       // Avoid division by -Inf
        p_no_outliers = (std::min) (1.0 - std::numeric_limits<double>::epsilon (), p_no_outliers);   // Avoid division by 0.
        k = log_probability / log(p_no_outliers);
      //std::cout<<"log_probability: "<<log_probability<<" log p_no_outliers: "<<p_no_outliers<<" k: "<<k<<" w: "<<w<<"inlier: "<<maybe_inliers_num<<std::endl;
      }
      else
      {
        delete [] cur_inti_sub_set_indexs;
        //--iterations;  
      }

      ++iterations;
    }

  std::set<int *, SubSetIndexComparator >::iterator it = chosen_sub_sets.begin();
  std::set<int *, SubSetIndexComparator >::iterator chosenSubSetsEnd = chosen_sub_sets.end();
  while(it!=chosenSubSetsEnd) {
    delete [] (*it);
    it++;
  }
  chosen_sub_sets.clear();


  std::vector<double> choosex(best_model_num, 0);
  std::vector<double> choosey(best_model_num, 0);
  int k1=0;
  int k2=0;
  for(int j=0; j<num_data; j++) {
    if(best_votes[j]){
      choose[j]=1;
       //std::cout<<j<<" ";
      choosex[k1++]=x[j];
      choosey[k2++]=y[j];
      }
      else{
      choose[j]=0;
      }
  }
  mylr.caltheta(choosex,choosey,best_model_num);
  //std::cout<<"a: "<<mylr.theta1<<" b: "<<mylr.theta0<<std::endl;

  delete []best_votes;
  delete []cur_votes;


  return best_model_num;
}
