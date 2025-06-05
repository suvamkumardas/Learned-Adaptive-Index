#ifndef LINE_SEGMENTS_IMPL_HPP
#define LINE_SEGMENTS_IMPL_HPP

#include "BaseLearnedModel.hpp"
#include <vector>
#include <cmath>
#include <tuple>

template <typename KeyType> class listOfSegments;  

template <typename KeyType>
class segment
{
    private:
        KeyType startKey, endkey;
        long startPos, endPos;
        double slope;
        int error;

    public:
        segment(KeyType startKey, KeyType endkey, long startPos, long endPos, double slope, int error)
        {
            this->startKey = startKey;
            this->endkey = endkey;
            this->startPos = startPos;
            this->endPos = endPos;
            this->slope = slope;
            this->error = error;
        }

        friend class listOfSegments<KeyType>;    
};

template <typename KeyType>
class listOfSegments: public baseLearnedModel<KeyType>
{
    private:
        std::vector<KeyType> startKeyList, endKeyList;
        std::vector<segment<KeyType>> segmentList;
        typename std::vector<KeyType>::iterator arrayBeginIter;
        long errorBound;
    
    public:

        listOfSegments(long errorBound)
        {
            this->errorBound = errorBound;
        }

        void buildIndex(std::vector<KeyType> &data, long startPos, long endPos)
        {
            this->arrayBeginIter = data.begin();

            while(startPos < endPos)
            {
                std::tuple<segment<KeyType>, bool, long> tup = this->buildSingleSegment(data, startPos, endPos);
                segment<KeyType> s = std::get<0>(tup);
                bool flag = std::get<1>(tup);
                long pos = std::get<2>(tup);

                segmentList.push_back(s);
                if (flag)
                    break;
                else
                    startPos = pos;    
            }

            if (startPos == endPos)
            {
                segment<KeyType> s = segmentList.back();
                KeyType lastKey = data[startPos];
                long startPos = s.startPos;
                KeyType startKey = s.startKey;
                double finalSlope = ((double)(endPos - startPos))/ ((double)(lastKey - startKey));
                long maxError = this->getErrorBound(data, startPos, startKey, endPos, finalSlope);
                s.slope = finalSlope;
                s.error = maxError;
                s.endPos = endPos;
                s.endkey = lastKey;

                //segmentList.push_back(s);                
            }

            startKeyList.reserve(segmentList.size());
            startKeyList.clear();
            endKeyList.reserve(segmentList.size());
            endKeyList.clear();

            for (segment<KeyType> s : segmentList)
            {
                startKeyList.push_back(s.startKey);
                endKeyList.push_back(s.endkey);
            }

        }

        long getPosition(KeyType key)
        {
            long pos = this->getSegmentPosForKey(key);
            segment<KeyType> s = segmentList[pos];
            double slope = s.slope;
            KeyType startKey = s.startKey, endKey = s.endkey;
            long startPos = s.startPos, endPos = s.endPos;
            int error = s.error;

            if (key == startKey)
                pos = startPos;
            else if (key == endKey)
                pos = endPos;
            else
            {
                long approxPos = round(slope * (key - startKey) + startPos);
                long approxPosMinusError =  approxPos - error;
                long approxPosPlusError  =  approxPos + error;
                if (approxPosMinusError < startPos)
                    approxPosMinusError = startPos;
                if (approxPosPlusError > endPos)
                    approxPosPlusError = endPos;                
                auto iter = std::lower_bound(this->arrayBeginIter + approxPosMinusError, this->arrayBeginIter + approxPosPlusError+1, key);
                pos = iter - this->arrayBeginIter;
            }        

            return pos;
        }

        long getSegmentPosForKey(KeyType key)
        {
            std::vector<long> begin(this->startKeyList.begin(), this->startKeyList.end());
            std::vector<long> end(this->endKeyList.begin(), this->endKeyList.end());
            std::vector<bool> temp;

            for(auto& element : begin)
                element = key - element;

            for(auto& element : end)
                element = element - key;

            temp.reserve(begin.size());
            temp.clear();
            for (int i = 0 ; i < begin.size() ; i++)
            {
                bool c = std::signbit(begin[i]) != std::signbit(end[i]);
                temp.push_back(c);
            }

            long pos;
            if (std::all_of(temp.cbegin(), temp.cend(), [](bool i) { return i ; }))
            {
                std::cout << "Error in finding segments" << std::endl;
                exit(0); 
            }
            else
            {
                for (int i = 0 ; i < temp.size() ; i++)
                    if (!temp[i]) {
                        pos = i;
                        break;
                    }
            }

            return pos;

        }

        std::tuple<segment<KeyType>, bool, long> buildSingleSegment(std::vector<KeyType> &data, long startPos, long endPos)
        {
            bool flag = true;
            KeyType startKey = data[startPos];
            double lowSlope = 0, highSlope = INFINITY;

            long pos;
            for(pos = startPos + 1 ; pos <= endPos ; pos++)
            {
                KeyType key = data[pos];
                double slope = ((double)(pos - startPos))/((double)(key - startKey));
                if (lowSlope <= slope && slope <= highSlope)
                {
                    long posPlusError = pos + this->errorBound;
                    long posMinusError = pos - this->errorBound;

                    if (posPlusError > endPos)
                        posPlusError = endPos;
                    if (posMinusError < startPos)
                        posMinusError = startPos;
                                            
                    double newLowSlope = ((double)(posMinusError - startPos))/((double)(key - startKey));
                    double newHighSlope = ((double)(posPlusError - startPos))/((double)(key - startKey));
                    if (newLowSlope > lowSlope)
                        lowSlope = newLowSlope;
                    if (newHighSlope < highSlope)
                        highSlope = newHighSlope;
                }
                else
                {
                    flag = false;
                    break;
                }       

            }

            //long lastPos = (!flag) ? pos - 1 : pos;
            long lastPos = pos - 1;
            KeyType lastKey = data[lastPos];
            double finalSlope = ((double)(lastPos - startPos))/ ((double)(lastKey - startKey));
            
            segment<KeyType> s(startKey, lastKey, startPos, lastPos, finalSlope, this->errorBound);
            std::tuple<segment<KeyType>, bool, long> tup{s, flag, pos};
            return tup;
        }

        long getErrorBound(std::vector<KeyType> &data, long startPos, KeyType startKey, long endPos, double slope)
        {
            long maxError = -1;
            long pos = startPos;
            while (pos <= endPos)
            {
                KeyType val = data[pos];
                long predPos = round(slope * (val - startKey) + startPos);
                long error = abs(pos - predPos);
                maxError = std::max(maxError, error);
                pos += 1;
            }
        
            return maxError;
        }
      
};

#endif