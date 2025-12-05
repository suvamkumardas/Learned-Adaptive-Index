#ifndef RADIXSPLINEIMPL_HPP
#define RADIXSPLINEIMPL_HPP

#include "../RadixSplineExtHash/include/rs/builder.h"
#include "BaseLearnedModel.hpp"
#include <vector>

template <typename KeyType>
class radixSplineModel: public baseLearnedModel<KeyType>
{
    private:
        std::vector<rs::Coord<KeyType>> &data;
        long startPos, endPos;
        rs::RadixSpline<KeyType> rs;
        //typename std::vector<rs::Coord<KeyType>>::iterator arrayBeginIter, globalStartIter;
        size_t initial_num_radix, max_error, bucket_size;

        static bool compare_by_key (const rs::Coord<KeyType>& p1, const rs::Coord<KeyType>& p2) {
            return p1.x < p2.x;
        }


    public:
        radixSplineModel(std::vector<rs::Coord<KeyType>> &inputData, size_t initial_num_radix = 18, 
            size_t error_bound = 32, size_t bucket_size = 1024) : data(inputData)
        {
            //this->data = data;
            this->initial_num_radix = initial_num_radix;
            this->max_error = error_bound;
            this->bucket_size = bucket_size;
        }
        
        void buildIndex(long startPos, long endPos)
        {
            KeyType min = data[startPos].x, max = data[endPos].x;
            rs::Builder<KeyType> rsb(min, max, this->initial_num_radix, this->max_error, this->bucket_size);
            for (int i = startPos ; i <= endPos ; i++) {
                const auto& key = data[i].x;
                rsb.AddKey(key);
            }
            this->rs = rsb.Finalize();
            //this->arrayBeginIter = std::begin(data);
            //this->globalStartIter = std::begin(data) + startPos;
            this->startPos = startPos;
            this->endPos = endPos;
        }

        void setStartAndEndPos(long startPos, long endPos)
        {
            this->startPos = startPos;
            this->endPos = endPos;
        }

        long getPositionGreaterThanOrEquals(KeyType key)
        {
            long pos = this->rs.getPositionGreaterThanOrEquals(data, this->startPos, key);
            return pos;
        }

        std::pair<rs::Coord<KeyType>*, size_t> rangeSearchOnInsertedData(KeyType low, KeyType high)
        {
            std::pair<rs::Coord<KeyType>*, size_t> res = this->rs.RangeSearch(data, low, high, true);
            return res;
        }

        int insertData(rs::Coord<KeyType> val)
        {
            return this->rs.insertNewKey(val.x, val.y);
        }

        std::pair<rs::Coord<KeyType>*, size_t> getAllElementsLessThanMin()
        {
            return this->rs.getAllElementsLessThanMin();
        }      

        std::pair<rs::Coord<KeyType>*, size_t> getAllElementsGreaterThanMax()
        {
            return this->rs.getAllElementsGreaterThanMax();
        }

        std::pair<long, long> retrainModel()
        {
            std::pair<rs::Coord<KeyType>*, size_t> p = this->rs.getAllInsertedData();
            rs::Coord<KeyType>* insertedData = p.first;
            size_t insertedDataSize = p.second;
            std::sort(insertedData, insertedData+insertedDataSize, compare_by_key);
            size_t newSize = this->data.size() + insertedDataSize;
            std::vector<rs::Coord<KeyType>> newData(newSize);
            std::copy(this->data.begin(), this->data.begin() + startPos, newData.begin());
            std::merge(this->data.begin() + startPos, this->data.begin() + endPos + 1, insertedData, insertedData+insertedDataSize, newData.begin() + startPos, compare_by_key);
            std::copy(this->data.begin() + endPos + 1, this->data.end(), newData.begin() + endPos + insertedDataSize + 1);        
            data = newData;
            buildIndex(startPos, endPos + insertedDataSize);
            return {startPos, endPos};
        }

        std::pair<long, long> getStartAndEndPos()
        {
            return {this->startPos, this->endPos};
        }      
};

#endif