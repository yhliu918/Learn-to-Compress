//============================================================================
// Name : Main.cpp
// Author : David Nogueira
//============================================================================

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <array>
#include <algorithm>
#include "headers/microunit.h"
#include "headers/easylogging++.h"
#include "headers/MLP.h"
INITIALIZE_EASYLOGGINGPP


// Example illustrating practical use of this MLP lib.
// Disclaimer: This is NOT an example of good machine learning practices
//              regarding training/testing dataset partitioning.

const int input_size = 6;
const int number_classes = 4;

const char *dataset = "out.txt";
const std::string mlp_weights = "classifier.mlp";


bool load_data(int *samples,
  std::vector<double>  *input,
  std::vector<double> *myclass) {
  // Load the block segment data-set. 
  std::ifstream srcFile("out.txt",std::ios::in); 
  if (!srcFile) {
    LOG(ERROR) << "Could not open file: " << dataset << ".";
    return false;
  }
  int lines=0;
  while(1){
    lines++;
    double tmp=0;
    int class_id=0;
    if(srcFile.eof()){break;}
    if(lines >=214300){break;}
    
    for (int j = 0; j < input_size; ++j) {
      srcFile>>tmp;
      
      input->push_back(tmp);
    }
     srcFile>>class_id; 
    for(int i=0;i<class_id;i++){
        myclass->push_back(0.0);
    }
    myclass->push_back(1.0);
    for(int j=class_id+1;j<number_classes;j++){
        myclass->push_back(0.0);
    }
    
  }
  LOG(INFO) << "Load "<<lines<<" into training set";
  (*samples)=lines;
  srcFile.close();
  return true;
}


int main(int argc, char *argv[]) {
  LOG(INFO) << "Test MLP with my dataset using backpropagation.";
  int samples = 0;
  std::vector<double> input;
  std::vector<double> myclass;

  // Load the data from file.
  if (!load_data(&samples, &input, &myclass)) {
    LOG(ERROR) << "Error processing input file.";
    return -1;
  }

  std::vector<TrainingSample> training_set;
  for (int j = 0; j < samples; ++j) {
    std::vector<double> training_set_input;
    std::vector<double> training_set_output;
    training_set_input.reserve(input_size);
    for (int i = 0; i < input_size; i++)
      training_set_input.push_back(*(&(input[0]) + j * input_size + i));
    training_set_output.reserve(number_classes);
    for (int i = 0; i < number_classes; i++)
      training_set_output.push_back(*(&(myclass[0]) + j * number_classes + i));
    training_set.emplace_back(std::move(training_set_input),
      std::move(training_set_output));
  }
  std::vector<TrainingSample> training_sample_set_with_bias(std::move(training_set));
  //set up bias
  for (auto & training_sample_with_bias : training_sample_set_with_bias) {
    training_sample_with_bias.AddBiasValue(1);
  }

  {
    // 4 inputs + 1 bias.
    // 1 hidden layer(s) of 4 neurons.
    // 3 outputs (1 per iris_class)
    MLP my_mlp({ input_size + 1, 50 ,number_classes }, { "sigmoid", "linear" }, false);

    int loops = 5000;

    // Train the network with backpropagation.
    LOG(INFO) << "Training for " << loops << " loops over data.";
    my_mlp.Train(training_sample_set_with_bias, .01, loops, 0.10, false);

    my_mlp.SaveMLPNetwork(mlp_weights);
  }
  
  //Destruction/Construction of a MLP object to show off saving and loading a trained model
    
  {
    MLP my_mlp(mlp_weights);

    int correct = 0;
    for (int j = 0; j < samples; ++j) {
      std::vector<double> guess;
      my_mlp.GetOutput(training_sample_set_with_bias[j].input_vector(), &guess);
      size_t class_id;
      my_mlp.GetOutputClass(guess, &class_id);

      if (myclass[j * number_classes + 0] == 1.0 && class_id == 0) {
        ++correct;
      }
      else if (myclass[j * number_classes + 1] == 1.0  && class_id == 1) {
        ++correct;
      }
      else if (myclass[j * number_classes + 2] == 1.0 && class_id == 2) {
        ++correct;
      }
     else if (myclass[j * number_classes + 3] == 1.0 && class_id == 3) {
        ++correct;
      }
    }
    LOG(INFO) << correct << "/" << samples
      << " (" << ((double)correct / samples * 100.0) << "%).";
  }
  return 0;
}
