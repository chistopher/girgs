
#include <iostream>

#include <Generator.h>


using namespace std;

int main(int argc, char* argv[])
{
    auto g = Generator(3);

    cout << sizeof(g) << endl;
    return 0;
}
