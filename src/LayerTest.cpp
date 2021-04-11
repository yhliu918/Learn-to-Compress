//============================================================================
// Name : LayerTest.cpp
// Author : David Nogueira
//============================================================================

#include "../headers/Layer.h"
#include "../headers/Sample.h"
#include "../headers/MLPUtils.h"

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>
#include "../headers/microunit.h"
#include "../headers/easylogging++.h"

INITIALIZE_EASYLOGGINGPP

int main(int argc, char* argv[]) {
  START_EASYLOGGINGPP(argc, argv);
  microunit::UnitTester::Run();
  return 0;
}