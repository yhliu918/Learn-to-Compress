#include <iostream>
#include "../headers/LinearRegression.h"
#include "../headers/Utils.h"

using namespace std;


void LinearRegression::train(double x[], double y[], int m ,double alpha, int iterations) {
    double *J = new double[iterations];
    this->theta = gradient_descent(x, y, alpha, iterations, J, m);
    /*
    cout << "J = ";
    for (int i = 0; i < iterations; ++i) {
        cout << J[i] << ' ';
    }
    cout << endl << "Theta: " << theta[0] << " " << theta[1] << endl;
    */
}

double LinearRegression::predict(double x) {
    return h(x, theta);
}

double LinearRegression::compute_cost(double x[], double y[], double theta[], int m) {
    double *predictions = calculate_predictions(x, theta, m);
    double *diff = Utils::array_diff(predictions, y, m);
    double *sq_errors = Utils::array_pow(diff, m, 2);
    return (1.0 / (2 * m)) * Utils::array_sum(sq_errors, m);
}

double LinearRegression::h(double x, double theta[]) {
    return theta[0] + theta[1] * x;
}

double *LinearRegression::calculate_predictions(double x[], double theta[], int m) {
    double *predictions = new double[m];

    // calculate h for each training example
    for (int i = 0; i < m; ++i) {
        predictions[i] = h(x[i], theta);
    }

    return predictions;
}

double *LinearRegression::gradient_descent(double x[], double y[], double alpha, int iters, double *J, int m) {
    double *theta = new double[2];
    theta[1] = (double) (y[m-1]-y[0])/ (x[m-1]-x[0]);
    theta[0] = y[0]-theta[1]*x[0];
    //cout << endl << "Theta: " << theta[0] << " " << theta[1] << endl;
    for (int i = 0; i < iters; ++i) {
        double *predictions = calculate_predictions(x, theta, m);
        double *diff = Utils::array_diff(predictions,y, m);
        /*
        for(int k=0;k<m;k++){
        cout<<diff[k]<<" ";
        }
        cout<<endl;
        */
        double *errors_x1 = diff;
        double *errors_x2 = Utils::array_multiplication(diff, x, m);

        theta[0] = theta[0] - alpha * (1.0 / m) * Utils::array_sum(errors_x1, m);
        theta[1] = theta[1] - alpha * (1.0 / m) * Utils::array_sum(errors_x2, m);
        //cout<< "new theta: "<<  theta[0] << " " << theta[1] << endl;
        //cout<<"errors :"<<alpha * (1.0 / m) * Utils::array_sum(errors_x1, m)<<" "<<alpha * (1.0 / m) * Utils::array_sum(errors_x2, m)<<endl;
        J[i] = compute_cost(x, y, theta, m);
    }

    return theta;
}

