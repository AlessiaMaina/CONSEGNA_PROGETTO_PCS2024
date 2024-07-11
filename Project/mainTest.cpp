#include <iostream>
#include "Tests.hpp"
#include "Utils.hpp"
#include "FracturesAndTraces.hpp"

using namespace std;

int main(int argc, char ** argv){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
