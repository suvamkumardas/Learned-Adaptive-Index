#include <iostream>
using namespace std;


enum queryType {
    BUILD_LEARNED_INDEX,
    QUERY_FROM_SAME_BOUND,
    QUERY_FROM_DIFFERENT_BOUND,
    CRACK_HIGH_VALUE,
    CRACK_LOW_VALUE,
    OVERLAP,
    MAX_QUERY_TYPE
};