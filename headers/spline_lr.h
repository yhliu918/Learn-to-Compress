
#ifndef SPLINE_LR_H_
#define SPLINE_LR_H_

#include "common.h"
#include "Utils.h"
#include "matrix_inverse.h"
#include<math.h>

//(A^TA)^(-1)A^TY
struct spline_lr{//theta + theta0*logx+ theta1*x + theta2* x^2 + theta3 * e^x
    double alpha;
    double theta1;
    double theta2;
    double theta3;

    
    
void caltheta(double x[], double y[], int m){
    double **matrix=new double*[m];
    for(int i=0;i<m;i++){
        matrix[i]= new double[4];
    }
    for(int i=0;i<m;i++){
        matrix[i][0]=1;
        /*
        if(x[i]>0.01){
            matrix[i][1]=log2(x[i]);
        }
        else{
            matrix[i][1]=0.0;
        }
        */
        matrix[i][1]=x[i];
        matrix[i][2]=x[i]*x[i];
        matrix[i][3]=x[i]*x[i]*x[i];
    }
    
    double **result=new double*[4];
    for(int i=0;i<4;i++){
        result[i]= new double[1];
    }

    
    double **ATA_reverse2=new double*[4];
    for(int i=0;i<4;i++){
        ATA_reverse2[i]= new double[4];
    }
    
    double **ATA=new double*[4];
    for(int i=0;i<4;i++){
        ATA[i]= new double[4];
    }
    
    double ATA2 [16]={0};

    double **AT=new double*[4];
    for(int i=0;i<4;i++){
        AT[i]= new double[m];
    }
    
    double **ATA_reverseAT=new double*[4];
    for(int i=0;i<4;i++){
        ATA_reverseAT[i]= new double[m];
    }

    AT=mat_T(4,m,matrix);
    ATA=multi(4,m,4,AT,matrix);
    
    
    for(int i=0;i<4;i++){
        for(int j=0;j<4;j++){
            ATA2[i*4+j] = ATA[i][j];
        }
    }
    
    /*
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            
            std::cout<<ATA[i][j]<<" ";
        }
        std::cout<<std::endl;
    }
    */
    double *ATA_reverse =new double[16];

    ATA_reverse=LUP_solve_inverse(ATA2);

    
    for(int i=0;i<4;i++){
        for(int j=0;j<4;j++){
            ATA_reverse2[i][j] = ATA_reverse[i*4+j];
            //std::cout<<ATA_reverse2[i][j]<<" ";
        }
        //std::cout<<std::endl;
    }
    ATA_reverseAT=multi(4,4,m,ATA_reverse2,AT);
    

    double **tmpy=new double*[m];
    for(int i=0;i<m;i++){
        tmpy[i]= new double[1];
        tmpy[i][0]=y[i];
    }
    

    result=multi(4,m,1,ATA_reverseAT,tmpy);

    alpha = result [0][0];
    theta1 = result [1][0];
    theta2 = result [2][0];
    theta3 = result [3][0];
    //std::cout<<"alpha: "<<alpha<<" theta1: "<<theta1<<" theta2: "<<theta2<<" theta3: "<<theta3<<std::endl;
    
}

    
};

    
/*

struct pol_lr{//theta + theta0*logx+ theta1*x + theta2* x^2 + theta3 * e^x
    double alpha;
    double theta1;
    double theta2;

    
    
void caltheta(double x[], double y[], int m){
    double **matrix=new double*[m];
    for(int i=0;i<m;i++){
        matrix[i]= new double[3];
    }
    for(int i=0;i<m;i++){
        matrix[i][0]=1;
        matrix[i][1]=x[i];
        matrix[i][2]=x[i]*x[i];
    }
    
    double **result=new double*[3];
    for(int i=0;i<3;i++){
        result[i]= new double[1];
    }

    
    double **ATA_reverse2=new double*[3];
    for(int i=0;i<3;i++){
        ATA_reverse2[i]= new double[3];
    }
    
    double **ATA=new double*[3];
    for(int i=0;i<3;i++){
        ATA[i]= new double[3];
    }
    
    double ATA2 [9]={0};

    double **AT=new double*[3];
    for(int i=0;i<3;i++){
        AT[i]= new double[m];
    }
    
    double **ATA_reverseAT=new double*[3];
    for(int i=0;i<3;i++){
        ATA_reverseAT[i]= new double[m];
    }

    AT=mat_T(3,m,matrix);
    ATA=multi(3,m,3,AT,matrix);
    
    
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            ATA2[i*3+j] = ATA[i][j];
        }
    }
    
    
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            
            std::cout<<ATA[i][j]<<" ";
        }
        std::cout<<std::endl;
    }
    
    double *ATA_reverse =new double[9];

    ATA_reverse=LUP_solve_inverse(ATA2);

    
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            ATA_reverse2[i][j] = ATA_reverse[i*3+j];
            //std::cout<<ATA_reverse2[i][j]<<" ";
        }
        //std::cout<<std::endl;
    }
    ATA_reverseAT=multi(3,3,m,ATA_reverse2,AT);
    

    double **tmpy=new double*[m];
    for(int i=0;i<m;i++){
        tmpy[i]= new double[1];
        tmpy[i][0]=y[i];
    }
    

    result=multi(3,m,1,ATA_reverseAT,tmpy);

    alpha = result [0][0];
    theta1 = result [1][0];
    theta2 = result [2][0];
    //std::cout<<"alpha: "<<alpha<<" theta1: "<<theta1<<" theta2: "<<theta2<<std::endl;
    
}

    
};

*/

#endif 
