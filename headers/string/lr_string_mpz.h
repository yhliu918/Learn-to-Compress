
#ifndef LR_STRING_MPZ_H_
#define LR_STRING_MPZ_H_
#include "common.h"
#include "Utils.h"
#include <math.h>

struct string_mpz_lr
{

    mpz_t theta0;
    mpz_t theta1;

    void caltheta(std::vector<int> &x, std::vector<mpz_t *> &y, int m)
    {

        mpz_init(theta0);
        mpz_init(theta1);

        mpz_t sumx;
        mpz_t sumy;
        mpz_t sumxy;
        mpz_t sumx2;
        mpz_init(sumx);
        mpz_init(sumy);
        mpz_init(sumxy);
        mpz_init(sumx2);

        for (int i = 0; i < m; i++)
        {
            //std::cout<<*(y[i])<<std::endl;
            mpz_add_ui(sumx, sumx, x[i]);
            mpz_add(sumy, sumy, *y[i]);
            mpz_addmul_ui(sumxy, *y[i], x[i]);
            mpz_add_ui(sumx2, sumx2, x[i] * x[i]);
        }

        // mpz_submul(avxy, avx, avy); // avxy = avxy - avx * avy = sumxy / m - avx * avy = (sumxy - avx * avy * m) / m = (sumxy - sumx / m * sumy / m * m) / m = (sumxy-sumx*sumy/m)/m= (sumxy*m-sumx*sumy)/m/m
        // (sumxy*m-sumx*sumy)/m/m
        
        mpz_t ccc;
        mpz_init(ccc);
        mpz_mul_ui(ccc, sumxy, m);
        mpz_submul(ccc, sumx, sumy);


        mpz_t xxx;
        mpz_init(xxx);
        mpz_mul_ui(xxx, sumx2, m);
        mpz_submul(xxx, sumx, sumx);
        

        // mpz_submul(avx2, avx, avx); // avx2 = avx2 - avx * avx = sumx2/m - sumx*sumx/m/m = (sumx2*m-sumx*sumx)/m/m)
        // (sumx2*m-sumx*sumx)/m/m)

        mpz_div(theta1, ccc, xxx);
        // theta0 * m = sumy - theta1 * sumx;
        mpz_mul(theta0, theta1, sumx);
        mpz_sub(theta0, sumy, theta0);
        mpz_div_ui(theta0, theta0, m);

        mpz_clear(sumx);
        mpz_clear(sumy);
        mpz_clear(sumxy);
        mpz_clear(sumx2);
        mpz_clear(ccc);
        mpz_clear(xxx);
    }
};

#endif
