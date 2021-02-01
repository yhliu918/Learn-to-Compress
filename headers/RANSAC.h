
#ifndef RANSAC_H_
#define RANSAC_H_

#include <set>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include<iostream>
#include "lr.h"






class RANSAC 
{
public:
  RANSAC();
  int compute(lr& mylr, double x[], double y[],int numbers,int num_for_estimate,bool *choose);


private:

  class SubSetIndexComparator 
  {
  private:
    int m_length;
  public:
    SubSetIndexComparator(int arrayLength) : m_length(arrayLength){}
    bool operator()(const int *arr1, const int *arr2) const 
    {
      for(int i=0; i<m_length; i++)
      {  if(arr1[i] < arr2[i])
        {
          return true;
        }
        if (arr1[i] > arr2[i])
        {
          return false;
        }
      }
      return false;      
    }
  };

  double probability_;  
  int max_iterations_;
  
};


#endif