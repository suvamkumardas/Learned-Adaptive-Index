#ifndef RADIXSPLINEIMPL_HPP
#define RADIXSPLINEIMPL_HPP

#include "../RadixSpline/include/rs/builder.h"
#include "BaseLearnedModel.hpp"
#include <vector>

template <typename KeyType>
class radixSplineModel: public baseLearnedModel<KeyType>
{
    private:
        rs::RadixSpline<KeyType> rs;
        KeyType *arrayBeginIter, *globalStartIter;
        size_t max_error;

    public:
        radixSplineModel(size_t error_bound = 32)
        {
            this->max_error = error_bound;
        }
        
        void buildIndex(KeyType data[], long startPos, long endPos)
        {
            KeyType min = data[startPos], max = data[endPos];
            rs::Builder<KeyType> rsb(min, max, 18, this->max_error);
            for (int i = startPos ; i <= endPos ; i++) {
                const auto& key = data[i];
                rsb.AddKey(key);
            }
            this->rs = rsb.Finalize();
            this->arrayBeginIter = data;
            this->globalStartIter = data + startPos;
        }

        long getPosition(KeyType key)
        {
            rs::SearchBound bound = this->rs.GetSearchBound(key);
            auto start = this->globalStartIter + bound.begin, last = this->globalStartIter + bound.end + 1;
            long pos = std::lower_bound(start, last, key) - this->arrayBeginIter;
            return pos;
        }      
};

#endif