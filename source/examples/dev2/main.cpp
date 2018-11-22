
#include <iostream>

#include <girgs/girgs-version.h>


using namespace std;


int main(int argc, char* argv[]) {

    // Library name
    std::cout << "Library girgs::girgs" << std::endl;
    std::cout << "========================================" << std::endl;

    // Library version
    std::cout << "Version: " << GIRGS_VERSION << std::endl;
    std::cout << std::endl;

    // Library type (static or dynamic)
    #ifdef BASELIB_STATIC_DEFINE
        std::cout << "Library type: STATIC" << std::endl;
    #else
        std::cout << "Library type: SHARED" << std::endl;
    #endif

    return 0;
}
