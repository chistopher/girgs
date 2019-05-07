
#include <gmock/gmock.h>
#include <cmath>

double distance(const std::vector<double>& a, const std::vector<double>& b) {
    assert(a.size() == b.size());
    auto result = 0.0;
    for(auto d=0u; d<a.size(); ++d){
        auto dist = std::abs(a[d] - b[d]);
        dist = std::min(dist, 1.0-dist);
        result = std::max(result, dist);
    }
    return result;
}


int main(int argc, char* argv[])
{
    ::testing::InitGoogleMock(&argc, argv);
    return RUN_ALL_TESTS();
}
