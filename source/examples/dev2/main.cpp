
#include <iostream>
#include <vector>
#include <algorithm>

#include <girgs/girgs-version.h>


using namespace std;


int main(int argc, char* argv[]) {


    auto vec = vector<int>(100, 0);

    vec[51] = 99999;

    bool found  = find(vec.begin(), vec.end(), 99999) != vec.end();


    return 0;
}
