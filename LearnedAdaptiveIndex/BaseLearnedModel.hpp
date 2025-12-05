#ifndef BASE_LEARNED_MODEL
#define BASE_LEARNED_MODEL

//#include <vector>
#include "../RadixSplineExtHash/include/rs/common.h"

template <typename KeyType>
class baseLearnedModel
{
    public:
        virtual void buildIndex(long startPos, long endPos) = 0;
        virtual long getPositionGreaterThanOrEquals(KeyType key) = 0;
        virtual std::pair<rs::Coord<KeyType>*, size_t> rangeSearchOnInsertedData(KeyType low, KeyType high) = 0;
        virtual int insertData(rs::Coord<KeyType> val) = 0;
        virtual std::pair<rs::Coord<KeyType>*, size_t> getAllElementsLessThanMin() = 0;
        virtual std::pair<rs::Coord<KeyType>*, size_t> getAllElementsGreaterThanMax() = 0;
        virtual std::pair<long, long> retrainModel() = 0;
        virtual std::pair<long, long> getStartAndEndPos() = 0;
        virtual void setStartAndEndPos(long startPos, long endPos) = 0;
};
#endif