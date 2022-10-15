
#ifndef MATRIX_IN_sizeVERSE
#define MATRIX_IN_sizeVERSE
#include "common.h"

#define N_size 3
double** mat_T(int m,int n, double **a){
    double **result=new double*[m];
        for(int i=0;i<m;i++){
            result[i]= new double[n];
        }
    
    for(int i=0;i<m;i++){//m is the column number of a
        for(int j=0;j<n;j++){
            result[i][j]=a[j][i];
        }
    }
    return result;
}
double** multi(int m, int n,int q, double **a, double **b )
{
        double **result=new double*[m];
        for(int i=0;i<m;i++){
            result[i]= new double[q];
            
        }
    
        for (int i = 0; i < m; i++)  
        {
            for (int j = 0; j < q; j++)      
            {
                for (int k = 0; k < n; k++)  
                {

                    result[i][j] += a[i][k] * b[k][j];
                    
                }
            }
        }
    return result;
}








//LUP分解
void LUP_Descomposition(double A[N_size*N_size],double L[N_size*N_size],double U[N_size*N_size],int P[N_size])
{
    int row=0;
    for(int i=0;i<N_size;i++)
    {
        P[i]=i;
    }
    for(int i=0;i<N_size-1;i++)
    {
        double p=0.0;
        for(int j=i;j<N_size;j++)
        {
            if(abs(A[j*N_size+i])>p)
            {
                p=abs(A[j*N_size+i]);
                row=j;
            }
        }
        if(0==p)
        {
            std::cout<< "矩阵奇异，无法计算逆" <<std::endl;
            return ;
        }

        //交换P[i]和P[row]
        int tmp=P[i];
        P[i]=P[row];
        P[row]=tmp;

        double tmp2=0.0;
        for(int j=0;j<N_size;j++)
        {
            //交换A[i][j]和 A[row][j]
            tmp2=A[i*N_size+j];
            A[i*N_size+j]=A[row*N_size+j];
            A[row*N_size+j]=tmp2;
        }

        //以下同LU分解
        double u=A[i*N_size+i],l=0.0;
        for(int j=i+1;j<N_size;j++)
        {
            l=A[j*N_size+i]/u;
            A[j*N_size+i]=l;
            for(int k=i+1;k<N_size;k++)
            {
                A[j*N_size+k]=A[j*N_size+k]-A[i*N_size+k]*l;
            }
        }

    }

    //构造L和U
    for(int i=0;i<N_size;i++)
    {
        for(int j=0;j<=i;j++)
        {
            if(i!=j)
            {
                L[i*N_size+j]=A[i*N_size+j];
            }
            else
            {
                L[i*N_size+j]=1;
            }
        }
        for(int k=i;k<N_size;k++)
        {
            U[i*N_size+k]=A[i*N_size+k];
        }
    }

}

//LUP求解方程
double * LUP_Solve(double L[N_size*N_size],double U[N_size*N_size],int P[N_size],double b[N_size])
{
    double *x=new double[N_size]();
    double *y=new double[N_size]();

    //正向替换
    for(int i = 0;i < N_size;i++)
    {
        y[i] = b[P[i]];
        for(int j = 0;j < i;j++)
        {
            y[i] = y[i] - L[i*N_size+j]*y[j];
        }
    }
    //反向替换
    for(int i = N_size-1;i >= 0; i--)
    {
        x[i]=y[i];
        for(int j = N_size-1;j > i;j--)
        {
            x[i] = x[i] - U[i*N_size+j]*x[j];
        }
        x[i] /= U[i*N_size+i];
    }
    return x;
}

/*****************矩阵原地转置BEGIN_size********************/

int getN_sizeext(int i, int m, int n)
{
  return (i%n)*m + i/n;
}


int getPre(int i, int m, int n)
{
  return (i%m)*n + i/m;
}

/* 处理以下标i为起点的环 */
void movedata(double *mtx, int i, int m, int n)
{
  double temp = mtx[i]; // 暂存
  int cur = i;    // 当前下标
  int pre = getPre(cur, m, n);
  while(pre != i)
  {
    mtx[cur] = mtx[pre];
    cur = pre;
    pre = getPre(cur, m, n);
  }
  mtx[cur] = temp;
}

/* 转置，即循环处理所有环 */
void transpose(double *mtx, int m, int n)
{
  for(int i=0; i<m*n; ++i)
  {
    int next = getN_sizeext(i, m, n);
    while(next > i) 
      next = getN_sizeext(next, m, n);
    if(next == i)  
      movedata(mtx, i, m, n);
  }
}
/*****************矩阵原地转置EN_sizeD********************/

//LUP求逆(将每列b求出的各列x进行组装)
double * LUP_solve_inverse(double A[N_size*N_size])
{
    //创建矩阵A的副本，注意不能直接用A计算，因为LUP分解算法已将其改变
    double *A_mirror = new double[N_size*N_size]();
    double *inv_A=new double[N_size*N_size]();//最终的逆矩阵（还需要转置）
    double *inv_A_each=new double[N_size]();//矩阵逆的各列
    //double *B    =new double[N_size*N_size]();
    double *b    =new double[N_size]();//b阵为B阵的列矩阵分量

    for(int i=0;i<N_size;i++)
    {
        double *L=new double[N_size*N_size]();
        double *U=new double[N_size*N_size]();
        int *P=new int[N_size]();

        //构造单位阵的每一列
        for(int i=0;i<N_size;i++)
        {
            b[i]=0;
        }
        b[i]=1;

        //每次都需要重新将A复制一份
        for(int i=0;i<N_size*N_size;i++)
        {
            A_mirror[i]=A[i];
        }

        LUP_Descomposition(A_mirror,L,U,P);

        inv_A_each=LUP_Solve (L,U,P,b);
        memcpy(inv_A+i*N_size,inv_A_each,N_size*sizeof(double));//将各列拼接起来
    }
    transpose(inv_A,N_size,N_size);//由于现在根据每列b算出的x按行存储，因此需转置

    return inv_A;
}





#endif
