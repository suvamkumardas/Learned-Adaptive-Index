#ifndef LAI_HPP
#define LAI_HPP

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <functional> 
#include <chrono>
#include "RadixSplineImpl.hpp"
#include "LineSegmentsImpl.hpp"
#include "queryType.hpp"
#include "BaseLearnedModel.hpp"
#include "../RadixSplineExtHash/include/rs/common.h"

#define EXT_HASH_IMPL


template <typename KeyType>
class LAI{

    private:
        std::vector<rs::Coord<KeyType>> data;
        std::vector<KeyType> lowValTable;
        std::vector<KeyType> highValTable;
        std::vector<long> lowPosTable;
        std::vector<long> highPosTable;
        std::vector<baseLearnedModel<KeyType>* > learnedIndexTable;
        std::vector<rs::Coord<KeyType>> pendingInserts;
        int error;
        queryType qType;
        bool useLearnedSort, hasInserts;
        unsigned long queryTypeCount[MAX_QUERY_TYPE];
        double queryRunTime[MAX_QUERY_TYPE];
        size_t initial_num_radix, max_error, bucket_size;
        //double partitionTime, sortTime, learnedIndexBuildTime;

        static bool compare_by_key (const rs::Coord<KeyType>& p1, const rs::Coord<KeyType>& p2) {
            return p1.x < p2.x;
        }

          
    public:

        LAI(std::vector<rs::Coord<KeyType>> data, size_t initial_num_radix=18, size_t max_error=32, size_t bucket_size=1024, bool useLearnedSort=true)
        {
            this->data = data;
            this->data.insert(this->data.begin(), rs::Coord<KeyType>{static_cast<KeyType>(-INFINITY), -1});
            this->data.push_back(rs::Coord<KeyType>{static_cast<KeyType>(INFINITY), -1});
            this->initial_num_radix = initial_num_radix;
            this->max_error = max_error;
            this->bucket_size = bucket_size;
            this->useLearnedSort = useLearnedSort;
            //this->partitionTime = this->sortTime = this->learnedIndexBuildTime = 0;

            for (int i = 0 ; i < MAX_QUERY_TYPE ; i++)
            {
                this->queryTypeCount[i] = 0;
                this->queryRunTime[i] = 0;
            }
            //pendingInserts = new std::vector<rs::Coord<KeyType>>();
            pendingInserts.reserve(this->data.size()/2);
            hasInserts = false;

        }

        void updatePendingInserts(std::vector<rs::Coord<KeyType>> &newData)
        {
            hasInserts = true;
            if(this->pendingInserts.size() == 0)
            {
                pendingInserts.reserve(newData.size());
                for (rs::Coord<KeyType> data:newData)
                    pendingInserts.push_back(data);
            }
            else
            {
                size_t finalSize = pendingInserts.size() + newData.size();
                std::vector<rs::Coord<KeyType>> tempData;
                tempData.resize(finalSize);
                std::merge(pendingInserts.begin(), pendingInserts.end(), newData.begin(), newData.end(), tempData.begin(), compare_by_key);
                pendingInserts = tempData;
            }
        }

        std::pair<rs::Coord<KeyType>*, size_t> getPendingData(KeyType low, KeyType high)
        {
            rs::Coord<KeyType>* res = nullptr;
            size_t size = 0;
            rs::Coord<KeyType> c = {low, -1};
            auto lowItr = std::lower_bound(pendingInserts.begin(), pendingInserts.end(), c, compare_by_key);
            if (lowItr != pendingInserts.end())
            {
                c = {high, -1};
                auto highItr = std::lower_bound(lowItr, pendingInserts.end(), c, compare_by_key);
                if (highItr == pendingInserts.end() || (*highItr).x > high)
                    highItr--;
                size = highItr - lowItr + 1;
                res = new rs::Coord<KeyType> [size];
                std::copy(lowItr, highItr+1, res);
                pendingInserts.erase(lowItr,highItr+1);   
            }
            return std::pair<rs::Coord<KeyType>*, size_t> {res, size};
        }

        std::pair<rs::Coord<KeyType>*, size_t> query(KeyType lowVal, KeyType highVal, bool measure_stats=true)
        {
            int size = this->lowValTable.size();
            long lowValPosInPartitionTable = -1, highValPosInPartitionTable = -1, lowValPosInArray, highValPosInArray;

            for (int i = 0 ; i < size ; i++)
            {
                if (this->lowValTable[i] <= lowVal && lowVal <= this->highValTable[i])
                {
                    lowValPosInPartitionTable = i;
                    break;
                }      
            }

            for (int i = 0 ; i < size ; i++)
            {
                if (this->lowValTable[i] <= highVal && highVal <= this->highValTable[i])
                {
                    highValPosInPartitionTable = i;
                    break;
                }      
            }

            auto started = std::chrono::high_resolution_clock::now();  
            std::tuple<long, long, rs::Coord<KeyType>*, size_t> pos;
            if (lowValPosInPartitionTable == -1 && highValPosInPartitionTable == -1)
            {
                if (checkForOverlapQuery(lowVal, highVal))
                {
                    this->qType = OVERLAP;
                    pos = this->executeOverlapQuery(lowVal, highVal);
                }
                else
                {    
                    this->qType = BUILD_LEARNED_INDEX;
                    pos = this->buildLearnedIndex(lowVal, highVal);
                }    
            }
            else if (lowValPosInPartitionTable == highValPosInPartitionTable)
            {
                this->qType = QUERY_FROM_SAME_BOUND;
                pos = this->getResultsFromSameBound(lowVal, highVal, lowValPosInPartitionTable);
            }
            else if(lowValPosInPartitionTable >= 0 && highValPosInPartitionTable >= 0)
            {
                this->qType = QUERY_FROM_DIFFERENT_BOUND;
                pos = this->getResultsFromDifferentBounds(lowVal, highVal, lowValPosInPartitionTable, highValPosInPartitionTable);
            }
            else if(lowValPosInPartitionTable >= 0 && highValPosInPartitionTable == -1)
            {
                this->qType = CRACK_HIGH_VALUE;
                pos = this->crackForHighValue(lowVal, lowValPosInPartitionTable, highVal);
            }
            else
            {
                this->qType = CRACK_LOW_VALUE;
                pos = this->crackForLowValue(lowVal, highVal, highValPosInPartitionTable);
            }
            auto done = std::chrono::high_resolution_clock::now();
            double time = (double)std::chrono::duration_cast<std::chrono::nanoseconds>(done-started).count();

            if (measure_stats) {
                this->queryTypeCount[this->qType]++;
                this->queryRunTime[this->qType] += time;
                //std::cout << "\nTime = " << time;
                //std::cout.flush();
            }

            #ifdef DEBUG
            bool status = this->checkAnswers(lowVal, highVal, pos.first, pos.second);
            std::cout << "\nQuery = (" << lowVal << " , " << highVal << ") : ";
            std::cout << pos.first << " " << pos.second << endl;
            std::cout << status << endl;
            std::cout << "****************************************";
            #endif
            
            long lowPosInArray = std::get<0>(pos), highPosInArray = std::get<1>(pos);
            rs::Coord<KeyType> *resInBuffer = std::get<2>(pos);
            size_t sizeInBuffer = std::get<3>(pos);
            print_buffer(resInBuffer, sizeInBuffer);

            size_t sizeInArr = highPosInArray - lowPosInArray + 1;
            size_t totalSize = sizeInArr + sizeInBuffer;
            rs::Coord<KeyType> *finalResult = new rs::Coord<KeyType> [sizeInArr + sizeInBuffer];
            std::copy(this->data.begin() + lowPosInArray, this->data.begin() + highPosInArray + 1, finalResult);
            std::copy(resInBuffer, resInBuffer + sizeInBuffer, finalResult+sizeInArr);
            print_buffer(finalResult, totalSize);
            return {finalResult, totalSize};
            
        }


        
        long partition(long startPos, long endPos, KeyType pivot)
        {    
            long long pivotPos;

            //find whether the pivot is present in the sub-array or not
            auto startIter(data.begin() + startPos);
            auto endIter(data.begin() + endPos+1);
            //auto lower = std::find(startIter, endIter, pivot);
            //const bool found = lower != endIter && *lower == pivot; // check that value has been found
            const bool found = true;
            
            #if 0
            if(found)
            {
                auto idx = lower - startIter;
                pivotPos = startPos + idx;
            }
            #endif    
            
            long i = startPos - 1;
            long j = endPos + 1;
            
            while(i <= j)
            {
                while(true)
                {
                    i+= 1;
                    if (data[i].x > pivot)
                        break;
                    if (data[i].x == pivot)
                        pivotPos = i;    
                }                        
                while(true)
                {
                    j-= 1;
                    if (data[j].x <= pivot)
                        break;
                }        

                if (data[j].x == pivot)
                    pivotPos = j;

                if (i < j)
                { 
                    std::swap(data[i].x, data[j].x);
                       
                    if (data[i].x == pivot)
                        pivotPos = i;
                    else if (data[j].x == pivot)
                        pivotPos = j; 
                       
                }   

            }
            if (found)
            {
                std::swap(data[j].x, data[pivotPos].x);
            }
                
            return j;
        }    
        



        std::pair<long, long>  crackAndSort(KeyType lowVal, KeyType highVal, long startPos, long endPos)
        {
            if (startPos == -1 and endPos == -1) 
            {
                std::pair<long, long> p = getPositions(lowVal, highVal);
                startPos = p.first;
                endPos = p.second;
            }

            //auto started = std::chrono::high_resolution_clock::now();
            long p = this->partition(startPos, endPos, lowVal);
            long q = this->partition(p+1, endPos, highVal);
            //auto iter = std::partition(data.begin() + startPos, data.begin() + endPos + 1, [lowVal](KeyType k) {return k < lowVal;});
            //auto pos = std::find(data.begin() + startPos, data.begin() + endPos + 1, lowVal);
            //std::swap(*iter,*pos);
            //long p = iter - data.begin();

            //iter = std::partition(data.begin() + p+1, data.begin() + endPos + 1, [highVal](KeyType k) {return k < highVal;});
            //pos = std::find(data.begin() + p+1, data.begin() + endPos + 1, highVal);
            //std::swap(*iter,*pos);
            //long q = iter - data.begin();
            
            //auto done = std::chrono::high_resolution_clock::now();
            //this->partitionTime += (double)std::chrono::duration_cast<std::chrono::nanoseconds>(done-started).count();
            
            //auto pos = std::partition(this->data.begin()+startPos, this->data.begin()+endPos+1, [lowVal](KeyType x){return lowVal > x;});
            //auto pos = std::partition(this->data.begin()+startPos, this->data.begin()+endPos+1, std::bind2nd(less<KeyType>(), lowVal));
            //long p = pos - this->data.begin();
            //pos = std::partition(pos+1, this->data.begin()+endPos+1, [highVal](KeyType x){return highVal > x;});
            //pos = std::partition(this->data.begin()+startPos, this->data.begin()+endPos+1, std::bind2nd(less<KeyType>(), highVal));
            //long q = pos - this->data.begin();
            if (data[p].x != lowVal)
                p += 1;
        
            if (p < q)
            {
                auto startIter(data.begin() + p);
                auto endIter(data.begin() + q + 1);
                //started = std::chrono::high_resolution_clock::now(); 
                if(!this->useLearnedSort)
                    std::sort(startIter, endIter, compare_by_key);
                else
                    learnedSort(data, p, q);
                //done = std::chrono::high_resolution_clock::now();
                //this->sortTime += (double)std::chrono::duration_cast<std::chrono::nanoseconds>(done-started).count();    
            }
            return std::pair<long, long>(p,q);
        }
                

        /*
        long LAI :: getStartPos(long lowVal)
        {
            long pos = 1;
            for (const std::unordered_map<std::string, long> &d : segTable)
            {
                if (data[d.at("lowPos")] < lowVal)
                {
                    pos = d.at("highPos") + 1;
                    if (pos > data.size() - 2)
                        throw std::logic_error("Low value of Query out of array bounds");
                }        
                else
                    break;
            }
            return pos;    
        }


        long LAI :: getEndPos(long highVal)
        {
            long pos = data.size() - 2;
            for (long i = segTable.size() - 1; i >= 0; i--)
            {
                const std::unordered_map<std::string, long> &d = segTable[i];
                if (data[d.at("highPos")] > highVal)
                {
                    pos = d.at("lowPos") - 1;
                    if (pos < 1)
                        throw std::logic_error("High Value of Query out of array bound");
                }
                else
                    break;
            }
            return pos;
        }
        */            
            
        void learnedSort(std::vector<rs::Coord<KeyType>>&data, long startPos, long endPos)
        {
            #ifdef DEBUG_LEARNED_SORT
            vector<KeyType> data_2(data.begin(), data.end());
            auto tempIter1 = data_2.begin() + startPos + 1, tempIter2 = data_2.begin() + endPos;
            auto started = std::chrono::high_resolution_clock::now(); 
            std::sort(tempIter1, tempIter2);
            auto done = std::chrono::high_resolution_clock::now();
            double timeToNativeSort = (double)std::chrono::duration_cast<std::chrono::nanoseconds>(done-started).count();
            
            started = std::chrono::high_resolution_clock::now();
            #endif 
            
            rs::Coord<KeyType> startVal = data[startPos], endVal = data[endPos];
            double slope = (double)(endPos - startPos) / (double)(endVal.x - startVal.x);
            double intercept = endPos - (slope * endVal.x);

            vector<rs::Coord<KeyType>> overflow_buffer;
            long size = endPos - startPos - 1;
            overflow_buffer.reserve(size);
            overflow_buffer.clear();

            rs::Coord<KeyType> maxVal = data[data.size() - 1];
            vector<rs::Coord<KeyType>> temp(size, maxVal);

            //The elements at startPos and endPos are already in their correct locations, so no need
            //to consider them
            for (long i = startPos + 1 ; i < endPos ; i++)
            {
                rs::Coord<KeyType> key = data[i];
                long pos = round(slope * key.x + intercept) - (startPos + 1);
                if (temp[pos].x == maxVal.x)
                    temp[pos] = key;
                else
                {
                    overflow_buffer.push_back(key);
                }
            }

            auto new_end = std::remove(temp.begin(), temp.end(), maxVal);
            temp.erase(new_end, temp.end());
            //temp.erase(std::remove(temp.begin(), temp.end(), maxVal), temp.end());
            std::sort(overflow_buffer.begin(), overflow_buffer.end(), compare_by_key);
            vector<rs::Coord<KeyType>> temp2(size);
            
            auto startIter = data.begin() + startPos + 1;
            //std::merge(temp.begin(), temp.end(), overflow_buffer.begin(), overflow_buffer.end(), temp2.begin());

            //std::copy(temp2.begin(), temp2.end(), startIter);

            std::merge(temp.begin(), temp.end(), overflow_buffer.begin(), overflow_buffer.end(), startIter, compare_by_key);
            
            #ifdef DEBUG_LEARNED_SORT
            done = std::chrono::high_resolution_clock::now();
            double timeToLearnedSort = (double)std::chrono::duration_cast<std::chrono::nanoseconds>(done-started).count();

            std::cout << "\nNumber of Elements to sort = " << endPos - startPos + 1 << std::endl;
            std::cout << "Native Sort = " << timeToNativeSort << std::endl;
            std::cout << "Learned Sort = " << timeToLearnedSort << std::endl;
            #endif            
        }
        
        long  getStartPos(KeyType lowVal)
        {
            auto iter = std::upper_bound(this->lowValTable.begin(), this->lowValTable.end(), lowVal);
            long pos;
            if (iter == this->lowValTable.begin()) 
            {
                pos = 1;
                return pos;
            } 

            long partitionInfoPos;
            if (iter == this->lowValTable.end())
                partitionInfoPos = this->lowValTable.size() - 1;
            else
                partitionInfoPos = iter - this->lowValTable.begin() - 1;
        
            pos = this->highPosTable[partitionInfoPos] + 1;
            return pos;
        }


        
        long getEndPos(KeyType highVal)
        {
            auto iter = std::upper_bound(this->highValTable.begin(), this->highValTable.end(), highVal);
            long pos;
            if (iter == this->highValTable.end()) 
            {
                pos = this->data.size() - 2;
                return pos;
            } 

            long partitionInfoPos = iter - this->highValTable.begin();  
            pos = this->lowPosTable[partitionInfoPos] - 1;
            return pos;
        }            


        
        std::pair<long, long> getPositions(KeyType lowVal, KeyType highVal)
        {
            long startPos;
            long endPos;
            if (this->lowValTable.size() == 0) 
            {
                startPos = 1;
                endPos = data.size() - 2;
            }
            else
            {     
                startPos = getStartPos(lowVal);
                endPos = getEndPos(highVal);
            }
            return {startPos, endPos};       
        }


        
        KeyType  getMaxElement(long start, long end)
        {
            KeyType max_elem = this->data[start].x;
            for (int i = start+1 ; i<= end ; i++)
            {
                KeyType elem = this->data[i].x;
                max_elem = std::max(max_elem, this->data[i].x);
            }
            return max_elem;
        }


        
        KeyType  getMinElement(long start, long end)
        {
            KeyType min_elem = this->data[start].x;
            for (int i = start+1 ; i<= end ; i++)
            {
                KeyType elem = this->data[i].x;
                min_elem = std::min(min_elem, this->data[i].x);
            }
            return min_elem;
        }


        
        std::pair<KeyType, KeyType>  getMinMaxElement(long start, long end)
        {
            KeyType min_elem = this->data[start].x, max_elem = this->data[start].x;
            for (int i = start+1 ; i<= end ; i++)
            {
                KeyType elem = this->data[i].x;
                min_elem = std::min(min_elem, this->data[i].x);
                max_elem = std::max(max_elem, this->data[i].x);
            }
            return {min_elem, max_elem};
        }


        /*
        bool LAI :: checkForOverlapQuery(long lowVal, long highVal)
        { 
            for (int i = 0 ; i < this->lowValTable.size() ; i++)
                if (lowVal < this->lowValTable[i] && highVal > this->highValTable[i])
                    return true;
            return false;
        }
        */


        
        bool  checkForOverlapQuery(KeyType lowVal, KeyType highVal)
        { 
            auto posIter = std::upper_bound (this->lowValTable.begin(), 
                                                        this->lowValTable.end(), lowVal);
            if (posIter == this->lowValTable.end())
                return false;
            
            long pos = posIter - this->lowValTable.begin();
            if (highVal > this->highValTable[pos])
                return true;
            
            return false;
        }


        
        baseLearnedModel<KeyType>*  buildLearnedIndexModel(long startPos, long endPos)
        {
            baseLearnedModel<KeyType> *learnedIndexModel = new radixSplineModel<KeyType>(this->error);
            //baseLearnedModel<KeyType> *learnedIndexModel = new listOfSegments<KeyType>(this->error);
            learnedIndexModel->buildIndex(this->data, startPos, endPos);
            return learnedIndexModel;
        } 


        
#ifndef EXT_HASH_IMPL
        std::pair<long, long>  buildLearnedIndex(KeyType lowVal, KeyType highVal, long lowValPos=-1, long highValPos=-1)
        {
            std::pair<long, long> pos = crackAndSort(lowVal, highVal, lowValPos, highValPos);
            long startPos = pos.first, endPos = pos.second;
            
            if (startPos > endPos)
            {
                cout << "No element in this range" << endl;
                return {-1, -1};
            }
                
            if ((startPos == endPos) && 
                    ((this->qType == CRACK_LOW_VALUE) || (this->qType == CRACK_HIGH_VALUE)
                    || (this->qType == BUILD_LEARNED_INDEX) || (this->qType == OVERLAP)))
                return {startPos, endPos};
                
            //auto started = std::chrono::high_resolution_clock::now(); 
            baseLearnedModel<KeyType> *learnedIndexModel = this->buildLearnedIndexModel(startPos, endPos);
            //auto done = std::chrono::high_resolution_clock::now();
            //this->learnedIndexBuildTime += (double)std::chrono::duration_cast<std::chrono::nanoseconds>(done-started).count();

            auto insertPosIter = std::upper_bound (this->lowValTable.begin(), this->lowValTable.end(), lowVal);
            this->lowValTable.insert(insertPosIter, lowVal);
            
            insertPosIter = std::lower_bound (this->lowValTable.begin(), this->lowValTable.end(), lowVal);
            long insertPos =  insertPosIter - this->lowValTable.begin();
            auto insertPosIter2 = this->highValTable.begin() + insertPos;
            this->highValTable.insert(insertPosIter2, highVal);

            auto insertPosIter3 = this->lowPosTable.begin() + insertPos;
            this->lowPosTable.insert(insertPosIter3, startPos);

            auto insertPosIter4 = this->highPosTable.begin() + insertPos;
            this->highPosTable.insert(insertPosIter4, endPos);

            auto insertPosIter5 = this->learnedIndexTable.begin() + insertPos;
            this->learnedIndexTable.insert(insertPosIter5, learnedIndexModel);

            return {startPos, endPos};
        }
#else

        std::tuple<long, long, rs::Coord<KeyType>*, size_t>  buildLearnedIndex(KeyType lowVal, KeyType highVal, long lowValPos=-1, long highValPos=-1)
        {
            std::pair<long, long> pos = crackAndSort(lowVal, highVal, lowValPos, highValPos);
            long startPos = pos.first, endPos = pos.second;
            
            if (startPos > endPos)
            {
                cout << "No element in this range" << endl;
                return {-1, -1, nullptr, 0};
            }
                
            if ((startPos == endPos) && 
                    ((this->qType == CRACK_LOW_VALUE) || (this->qType == CRACK_HIGH_VALUE)
                    || (this->qType == BUILD_LEARNED_INDEX) || (this->qType == OVERLAP)))
                return {startPos, endPos, nullptr, 0};
                
            //auto started = std::chrono::high_resolution_clock::now(); 
            baseLearnedModel<KeyType> *learnedIndexModel = this->buildLearnedIndexModel(startPos, endPos);
            std::pair<rs::Coord<KeyType>*, size_t> insertedAndPendingResults = getResultsFromInsertedAndPendingData(learnedIndexModel, lowVal, highVal);
            
            //auto done = std::chrono::high_resolution_clock::now();
            //this->learnedIndexBuildTime += (double)std::chrono::duration_cast<std::chrono::nanoseconds>(done-started).count();

            auto insertPosIter = std::upper_bound (this->lowValTable.begin(), this->lowValTable.end(), lowVal);
            this->lowValTable.insert(insertPosIter, lowVal);
            
            insertPosIter = std::lower_bound (this->lowValTable.begin(), this->lowValTable.end(), lowVal);
            long insertPos =  insertPosIter - this->lowValTable.begin();
            auto insertPosIter2 = this->highValTable.begin() + insertPos;
            this->highValTable.insert(insertPosIter2, highVal);

            auto insertPosIter3 = this->lowPosTable.begin() + insertPos;
            this->lowPosTable.insert(insertPosIter3, startPos);

            auto insertPosIter4 = this->highPosTable.begin() + insertPos;
            this->highPosTable.insert(insertPosIter4, endPos);

            auto insertPosIter5 = this->learnedIndexTable.begin() + insertPos;
            this->learnedIndexTable.insert(insertPosIter5, learnedIndexModel);

            return {startPos, endPos, insertedAndPendingResults.first, insertedAndPendingResults.second};
        }
#endif

        
        long  getFirstEntryGreaterThanLowVal(KeyType lowVal)
        {
            auto insertPosIter = std::upper_bound (this->lowValTable.begin(), this->lowValTable.end(), lowVal);
            long pos = insertPosIter -  this->lowValTable.begin();
            //return (pos > 0) ? pos : 0; 
            return pos;                                         
        }


        
        long  getLastEntryLessThanHighVal(KeyType highVal)
        {
            auto insertPosIter = std::upper_bound (this->lowValTable.begin(), this->lowValTable.end(), highVal);
            long pos = insertPosIter -  this->lowValTable.begin() - 1;
        /*
            if (pos < this->lowValTable.size() - 1)
                return pos;
            return this->lowValTable.size() - 1;  
        */
            return pos;                                        
        }

        void print_buffer(rs::Coord<KeyType>* buff, size_t size)
        {
            std::cout << "\n";
            for (int i = 0 ; i < size ; i++)
            {
                std::cout << buff[i].x << " ";
            }
            std::cout << "\n";
        }


        std::pair<rs::Coord<KeyType>*, size_t> getResultsFromInsertedAndPendingData(baseLearnedModel<KeyType> *learnedIndexModel, KeyType lowVal, KeyType highVal)
        {
            rs::Coord<KeyType>  *resInBuffer = nullptr;
            size_t resInBufferSize = 0;
            if (hasInserts)
            {
                std::pair<rs::Coord<KeyType>*, size_t> p1 = learnedIndexModel->rangeSearchOnInsertedData(this->data, lowVal, highVal);
                print_buffer(p1.first, p1.second);
                std::pair<rs::Coord<KeyType>*, size_t> p2 = getPendingData(lowVal, highVal);
                print_buffer(p2.first, p2.second);
                vector<rs::Coord<KeyType>> failedToInsert;
                failedToInsert.reserve(p2.second); 

                for(int i = 0 ; i < p2.second ; i++)
                {
                    rs::Coord<KeyType> c = p2.first[i];
                    int status = learnedIndexModel->insertData(c);
                    if (status == -1)
                        failedToInsert.push_back(c);
                }
                if (failedToInsert.size() > 0)
                    updatePendingInserts(failedToInsert);

                resInBufferSize = p1.second + p2.second;
                resInBuffer = new rs::Coord<KeyType> [resInBufferSize];
                std::copy(p1.first, p1.first + p1.second, resInBuffer);
                std::copy(p2.first, p2.first + p2.second, resInBuffer+p1.second);
                print_buffer(resInBuffer, resInBufferSize);
            }
            return {resInBuffer, resInBufferSize};

        }

        //insert elements to models falling between two sorted bounds
        std::pair<rs::Coord<KeyType>*, size_t> insertDataBetweenBoundaries(long lowValPosInPartitionTable, long highValPosInPartitionTable)
        {
            KeyType tempHigh = this->highValTable[lowValPosInPartitionTable];
            KeyType tempLow = this->lowValTable[highValPosInPartitionTable];
            vector<rs::Coord<KeyType>> failedToInsert;
            rs::Coord<KeyType>* res = nullptr;
            size_t resSize = 0;
            failedToInsert.reserve(100000); 
            while (tempHigh < tempLow)
            {
                lowValPosInPartitionTable++;
                KeyType nextLow = lowValTable[lowValPosInPartitionTable];
                std::pair<rs::Coord<KeyType>*, size_t> p = getPendingData(tempHigh, nextLow);

                rs::Coord<KeyType>* temp = new rs::Coord<KeyType> [resSize + p.second];
                std::copy(res, res+resSize, temp);
                std::copy(p.first, p.first+p.second,temp+resSize);
                resSize += p.second;
                res = temp;

                baseLearnedModel<KeyType> *learnedIndexModel1 = learnedIndexTable[lowValPosInPartitionTable - 1];
                baseLearnedModel<KeyType> *learnedIndexModel2 = learnedIndexTable[lowValPosInPartitionTable];
                for(int i = 0 ; i < p.second ; i++)
                {
                    rs::Coord<KeyType> c = p.first[i];
                    int status = learnedIndexModel1->insertData(c);
                    if (status == -1)
                    {
                        status = learnedIndexModel2->insertData(c);
                        if (status == -1)
                            failedToInsert.push_back(c);
                    }
                }
                tempHigh = highValTable[lowValPosInPartitionTable];
            }
            
            if (failedToInsert.size() > 0)
                updatePendingInserts(failedToInsert);

            return {res, resSize};    
        }


        std::tuple<long, long, rs::Coord<KeyType>*, size_t>  getResultsFromSameBound(KeyType lowVal, KeyType highVal, long posInPartitionTable)
        {
            baseLearnedModel<KeyType> *learnedIndexModel = this->learnedIndexTable[posInPartitionTable];
            long lowValPosInArray = learnedIndexModel->getPosition(this->data,lowVal);
            long highValPosInArray = learnedIndexModel->getPosition(this->data,highVal);
            std::pair<rs::Coord<KeyType>*, size_t> insertedAndPendingResults = getResultsFromInsertedAndPendingData(learnedIndexModel, lowVal, highVal);
            return {lowValPosInArray, highValPosInArray, insertedAndPendingResults.first, insertedAndPendingResults.second};
        }


#ifndef EXT_HASH_IMPL        
        std::pair<long, long>  getResultsFromDifferentBounds(KeyType lowVal, KeyType highVal, long lowValPosInPartitionTable, long highValPosInPartitionTable)
        {
            baseLearnedModel<KeyType> *learnedIndexModel = this->learnedIndexTable[lowValPosInPartitionTable];
            long lowValPosInArray = learnedIndexModel->getPosition(this->data,lowVal);
            
            learnedIndexModel = this->learnedIndexTable[highValPosInPartitionTable];
            long highValPosInArray = learnedIndexModel->getPosition(this->data, highVal);

            this->buildIndexForGap(lowValPosInPartitionTable, highValPosInPartitionTable);
            return {lowValPosInArray, highValPosInArray};
        }
#else        

        std::tuple<long, long, rs::Coord<KeyType>*, size_t>  getResultsFromDifferentBounds(KeyType lowVal, KeyType highVal, long lowValPosInPartitionTable, long highValPosInPartitionTable)
        {
            KeyType tempHigh = this->highValTable[lowValPosInPartitionTable];
            std::tuple<long, long, rs::Coord<KeyType>*, size_t> firstResults = getResultsFromSameBound(lowVal, tempHigh, lowValPosInPartitionTable);
            rs::Coord<KeyType>* firstInserted = std::get<2>(firstResults);
            size_t firstSize = std::get<3>(firstResults);

            KeyType tempLow = this->lowValTable[highValPosInPartitionTable];
            std::tuple<long, long, rs::Coord<KeyType>*, size_t> secondResults = getResultsFromSameBound(tempLow, highVal, highValPosInPartitionTable);
            rs::Coord<KeyType>* secondInserted = std::get<2>(secondResults);
            size_t secondSize = std::get<3>(secondResults);

            std::pair<rs::Coord<KeyType>*, size_t> middleResults = this->buildIndexForGap(lowValPosInPartitionTable, highValPosInPartitionTable);
            rs::Coord<KeyType>* middleInserted = middleResults.first;
            size_t middleSize = middleResults.second;

            //insert elements to models falling between two sorted bounds
            auto it = std::lower_bound(this->lowValTable.begin(), this->lowValTable.end(), tempLow);
            long newHighValPosInPartitionTable = it - this->lowValTable.begin();
            std::pair<rs::Coord<KeyType>*, size_t> insertedDataBetween = insertDataBetweenBoundaries(lowValPosInPartitionTable, newHighValPosInPartitionTable);
            rs::Coord<KeyType>* betweenData = insertedDataBetween.first;
            size_t betweenDataSize = insertedDataBetween.second;
            
            size_t finalInsertedAndPendingResultsSize = firstSize + secondSize + middleSize + betweenDataSize;
            rs::Coord<KeyType>* finalInsertedAndPendingResults = new rs::Coord<KeyType> [finalInsertedAndPendingResultsSize];
            std::copy(firstInserted, firstInserted+firstSize, finalInsertedAndPendingResults);
            std::copy(secondInserted, secondInserted+secondSize, finalInsertedAndPendingResults+firstSize);
            std::copy(middleInserted, middleInserted+middleSize, finalInsertedAndPendingResults+firstSize+secondSize);
            std::copy(betweenData, betweenData+betweenDataSize, finalInsertedAndPendingResults+firstSize+secondSize+middleSize);
            return {std::get<0>(firstResults), std::get<1>(secondResults), finalInsertedAndPendingResults, finalInsertedAndPendingResultsSize};
        }
#endif


#ifndef EXT_HASH_IMPL        
        std::pair<long, long>  crackForLowValue(KeyType lowVal, KeyType highVal, long highValPosInPartitionTable, bool computeHigh=true)
        {
            long highValPosInArray = -1;
            if (computeHigh)
            {
                baseLearnedModel<KeyType> *learnedIndexModel = this->learnedIndexTable[highValPosInPartitionTable];
                highValPosInArray = learnedIndexModel->getPosition(this->data, highVal);
            }

            long firstPosInPartitionTable = this->getFirstEntryGreaterThanLowVal(lowVal);
            this->buildIndexForGap(firstPosInPartitionTable, highValPosInPartitionTable);

            long startPos = (firstPosInPartitionTable == 0) ? 1 : this->highPosTable[firstPosInPartitionTable - 1] + 1;
            long endPos = this->lowPosTable[firstPosInPartitionTable] - 1;

            //long endVal = *std::max_element(this->data.begin() + startPos, this->data.begin() + endPos);
            KeyType endVal = this->getMaxElement(startPos, endPos);
            std::pair p = this->buildLearnedIndex(lowVal, endVal, startPos, endPos);
            long lowValPosInArray = p.first;

            //if one point is left to be cracked, add that point in the beginning and rebuild the whole learned model again
            if (lowVal == endVal)
            {
                //delete (learnedIndexModel);
                
                baseLearnedModel<KeyType> *learnedIndexModel = this->buildLearnedIndexModel(lowValPosInArray, this->highPosTable[firstPosInPartitionTable]);
                this->lowValTable[firstPosInPartitionTable] = lowVal;
                this->lowPosTable[firstPosInPartitionTable] = lowValPosInArray;
                this->learnedIndexTable[firstPosInPartitionTable] = learnedIndexModel;
            }

            return {lowValPosInArray, highValPosInArray};
        }
#else

        std::tuple<long, long, rs::Coord<KeyType>*, size_t>  crackForLowValue(KeyType lowVal, KeyType highVal, long highValPosInPartitionTable, bool computeHigh=true)
        {
            long highValPosInArray = -1;
        #if 0
            if (computeHigh)
            {
                baseLearnedModel<KeyType> *learnedIndexModel = this->learnedIndexTable[highValPosInPartitionTable];
                highValPosInArray = learnedIndexModel->getPosition(this->data, highVal);
            }
        #endif

            KeyType tempLow = this->lowValTable[highValPosInPartitionTable];
            std::tuple<long, long, rs::Coord<KeyType>*, size_t> secondResults = getResultsFromSameBound(tempLow, highVal, highValPosInPartitionTable);
            highValPosInArray = std::get<1>(secondResults);
            rs::Coord<KeyType>* secondInserted = std::get<2>(secondResults);
            size_t secondSize = std::get<3>(secondResults);


            long firstPosInPartitionTable = this->getFirstEntryGreaterThanLowVal(lowVal);
            std::pair<rs::Coord<KeyType>*, size_t> gapInsertedResult = this->buildIndexForGap(firstPosInPartitionTable, highValPosInPartitionTable);
            rs::Coord<KeyType>* middleInserted = gapInsertedResult.first;
            size_t middleSize = gapInsertedResult.second;

            long startPos = (firstPosInPartitionTable == 0) ? 1 : this->highPosTable[firstPosInPartitionTable - 1] + 1;
            long endPos = this->lowPosTable[firstPosInPartitionTable] - 1;

            //long endVal = *std::max_element(this->data.begin() + startPos, this->data.begin() + endPos);
            KeyType endVal = this->getMaxElement(startPos, endPos);
            std::tuple<long, long, rs::Coord<KeyType>*, size_t> firstResults;
            firstResults = this->buildLearnedIndex(lowVal, endVal, startPos, endPos);
           
            //if one point is left to be cracked, add that point in the beginning and rebuild the whole learned model again
            if (lowVal == endVal)
            {
                std::cout << "crackForLowValue - This case is not handled properly" << std::endl;
                std::cout << "Low = " << lowVal << " , High = " << highVal << std::endl;
                /*
                baseLearnedModel<KeyType> *learnedIndexModel = learnedIndexTable[firstPosInPartitionTable];
                KeyType tempLow2 = (firstPosInPartitionTable == 0) ? this->data.begin() ? lowValTable[firstPosInPartitionTable - 1];
                KeyType tempHigh2 = (firstPosInPartitionTable == HighValTable.size() - 1) ? this->data.back() ? HighValTable[firstPosInPartitionTable - 1];
                std::pair<rs::Coord<KeyType>*, size_t> dataToInsert = learnedIndexModel->rangeSearchOnInsertedData(this->data, tempLow2, tempHigh2);
                
                tempHigh2 = this->highValTable[firstPosInPartitionTable];
                long tempPos = this->highPosTable[firstPosInPartitionTable];
                auto it = lowValTable.begin() + firstPosInPartitionTable;
                lowValTable.erase(it);
                it = highValTable.begin() + firstPosInPartitionTable;
                highValTable.erase(it);
                it = lowPosTable.begin() + firstPosInPartitionTable;
                lowPosTable.erase(it);
                it = highPosTable.begin() + firstPosInPartitionTable;
                highPosTable.erase(it);
                it = learnedIndexTable.begin() + firstPosInPartitionTable;
                learnedIndexTable.erase(it);

                firstResults = this->buildLearnedIndex(lowVal, tempHigh2, startPos, tempPos);
                //need to add dataToInsert into the new model
                */
            }

            long lowValPosInArray = std::get<0>(firstResults);
            rs::Coord<KeyType>* firstInserted = std::get<2>(firstResults);
            size_t firstSize = std::get<3>(firstResults);


            //insert elements to models falling between two sorted bounds
            auto it = std::lower_bound(this->highValTable.begin(), this->highValTable.end(), highVal);
            long newHighValPosInPartitionTable = it - this->highValTable.begin();
            std::pair<rs::Coord<KeyType>*, size_t> insertedDataBetween = insertDataBetweenBoundaries(firstPosInPartitionTable, newHighValPosInPartitionTable);
            rs::Coord<KeyType>* betweenData = insertedDataBetween.first;
            size_t betweenDataSize = insertedDataBetween.second;

            size_t finalInsertedAndPendingResultsSize = firstSize + secondSize + middleSize + betweenDataSize;
            rs::Coord<KeyType>* finalInsertedAndPendingResults = new rs::Coord<KeyType> [finalInsertedAndPendingResultsSize];
            std::copy(firstInserted, firstInserted+firstSize, finalInsertedAndPendingResults);
            std::copy(secondInserted, secondInserted+secondSize, finalInsertedAndPendingResults+firstSize);
            std::copy(middleInserted, middleInserted+middleSize, finalInsertedAndPendingResults+firstSize+secondSize);
            std::copy(betweenData, betweenData+betweenDataSize, finalInsertedAndPendingResults+firstSize+secondSize+middleSize);

            return {lowValPosInArray, highValPosInArray, finalInsertedAndPendingResults, finalInsertedAndPendingResultsSize};
        }

#endif


        
#ifndef EXT_HASH_IMPL
        std::pair<long, long>  crackForHighValue(KeyType lowVal, long lowValPosInPartitionTable, KeyType highVal, bool computeLow=true)
        {
            long lowValPosInArray = -1;
            if (computeLow)
            {
                baseLearnedModel<KeyType> *learnedIndexModel = this->learnedIndexTable[lowValPosInPartitionTable];
                lowValPosInArray = learnedIndexModel->getPosition(this->data, lowVal);
            }

            long lastPosInPartitionTable = this->getLastEntryLessThanHighVal(highVal);
            long oldPartitionTableSize = this->lowValTable.size();
            this->buildIndexForGap(lowValPosInPartitionTable, lastPosInPartitionTable);
            long newPartitionTableSize = this->lowValTable.size();
            lastPosInPartitionTable += newPartitionTableSize - oldPartitionTableSize; 

            long startPos = this->highPosTable[lastPosInPartitionTable] + 1;
            long endPos = (lastPosInPartitionTable < this->lowPosTable.size() - 1) ? 
                            this->lowPosTable[lastPosInPartitionTable + 1] - 1 : this->data.size() - 2;

            KeyType startVal = this->getMinElement(startPos, endPos);
            //KeyType startVal = *std::min_element(data.begin() + startPos, data.begin() + endPos + 1);
            std:: pair p = this->buildLearnedIndex(startVal, highVal, startPos, endPos);
            long highValPosInArray = p.second;
            
            //if one point is left to be cracked, add that point to the existing partition info
            if (startVal == highVal)
            {
                //baseLearnedModel<KeyType> *learnedIndexModel = this->buildLearnedIndexModel(this->highPosTable[lastPosInPartitionTable], highValPosInArray);
                baseLearnedModel<KeyType> *learnedIndexModel = this->buildLearnedIndexModel(this->lowPosTable[lastPosInPartitionTable], highValPosInArray);
                this->highValTable[lastPosInPartitionTable] = highVal;
                this->highPosTable[lastPosInPartitionTable] = highValPosInArray;
                this->learnedIndexTable[lastPosInPartitionTable] = learnedIndexModel;
            }
                           
            return {lowValPosInArray, highValPosInArray};
        }
#else
        std::tuple<long, long, rs::Coord<KeyType>*, size_t>  crackForHighValue(KeyType lowVal, long lowValPosInPartitionTable, KeyType highVal, bool computeLow=true)
        {
            long lowValPosInArray = -1;
#if 0
            if (computeLow)
            {
                baseLearnedModel<KeyType> *learnedIndexModel = this->learnedIndexTable[lowValPosInPartitionTable];
                lowValPosInArray = learnedIndexModel->getPosition(this->data, lowVal);
            }
#endif            

            KeyType tempHigh = this->highValTable[lowValPosInPartitionTable];
            std::tuple<long, long, rs::Coord<KeyType>*, size_t> firstResults = getResultsFromSameBound(lowVal, tempHigh, lowValPosInPartitionTable);
            rs::Coord<KeyType>* firstInserted = std::get<2>(firstResults);
            size_t firstSize = std::get<3>(firstResults);
            lowValPosInArray = std::get<0>(firstResults);

            long lastPosInPartitionTable = this->getLastEntryLessThanHighVal(highVal);
            long oldPartitionTableSize = this->lowValTable.size();
            std::pair<rs::Coord<KeyType>*, size_t> middleResults = this->buildIndexForGap(lowValPosInPartitionTable, lastPosInPartitionTable);
            long newPartitionTableSize = this->lowValTable.size();
            lastPosInPartitionTable += newPartitionTableSize - oldPartitionTableSize; 
            rs::Coord<KeyType>* middleInserted = middleResults.first;
            size_t middleSize = middleResults.second;

            long startPos = this->highPosTable[lastPosInPartitionTable] + 1;
            long endPos = (lastPosInPartitionTable < this->lowPosTable.size() - 1) ? 
                            this->lowPosTable[lastPosInPartitionTable + 1] - 1 : this->data.size() - 2;

            KeyType startVal = this->getMinElement(startPos, endPos);
            //KeyType startVal = *std::min_element(data.begin() + startPos, data.begin() + endPos + 1);
            std::tuple<long, long, rs::Coord<KeyType>*, size_t> secondResults = this->buildLearnedIndex(startVal, highVal, startPos, endPos);
            long highValPosInArray = std::get<1>(secondResults);
            rs::Coord<KeyType>* secondInserted = std::get<2>(secondResults);
            size_t secondSize = std::get<3>(secondResults);

            lastPosInPartitionTable = std::lower_bound(highValTable.begin(), highValTable.end(), highVal) - highValTable.begin();
            std::pair<rs::Coord<KeyType>*, size_t> insertedDataBetween = insertDataBetweenBoundaries(lowValPosInPartitionTable, lastPosInPartitionTable);
            rs::Coord<KeyType>* betweenData = insertedDataBetween.first;
            size_t betweenDataSize = insertedDataBetween.second;

            //if one point is left to be cracked, add that point to the existing partition info
            if (startVal == highVal)
            {
                //baseLearnedModel<KeyType> *learnedIndexModel = this->buildLearnedIndexModel(this->highPosTable[lastPosInPartitionTable], highValPosInArray);
                std::cout << "crackForHighValue - This case is not handled properly" << std::endl;
                std::cout << "Low = " << lowVal << " , High = " << highVal << std::endl;
                /*
                baseLearnedModel<KeyType> *learnedIndexModel = this->buildLearnedIndexModel(this->lowPosTable[lastPosInPartitionTable], highValPosInArray);
                this->highValTable[lastPosInPartitionTable] = highVal;
                this->highPosTable[lastPosInPartitionTable] = highValPosInArray;
                this->learnedIndexTable[lastPosInPartitionTable] = learnedIndexModel;
                */
            }
            size_t finalInsertedAndPendingResultsSize = firstSize + secondSize + middleSize + betweenDataSize;
            rs::Coord<KeyType>* finalInsertedAndPendingResults = new rs::Coord<KeyType> [finalInsertedAndPendingResultsSize];
            std::copy(firstInserted, firstInserted+firstSize, finalInsertedAndPendingResults);
            std::copy(secondInserted, secondInserted+secondSize, finalInsertedAndPendingResults+firstSize);
            std::copy(middleInserted, middleInserted+middleSize, finalInsertedAndPendingResults+firstSize+secondSize);
            std::copy(betweenData, betweenData+betweenDataSize, finalInsertedAndPendingResults+firstSize+secondSize+middleSize);

                           
            return {lowValPosInArray, highValPosInArray, finalInsertedAndPendingResults, finalInsertedAndPendingResultsSize};
        }
#endif


#ifndef EXT_HASH_IMPL        
        std::pair<long, long>  executeOverlapQuery(KeyType lowVal, KeyType highVal)
        {
            long firstPosInPartitionTable = this->getFirstEntryGreaterThanLowVal(lowVal);
            long lastPosInPartitionTable = this->getLastEntryLessThanHighVal(highVal);
            KeyType tempHighVal = this->lowValTable[firstPosInPartitionTable];
            KeyType tempLowVal = this->highValTable[lastPosInPartitionTable];

            long oldPartitionTableSize = this->lowValTable.size();
            std::pair<long, long> pair = this->crackForLowValue(lowVal, tempHighVal, firstPosInPartitionTable, false);
            long lowValPosInArray = pair.first;

            long newPartitionTableSize = this->lowValTable.size();
            lastPosInPartitionTable += newPartitionTableSize - oldPartitionTableSize;
            this->buildIndexForGap(firstPosInPartitionTable, lastPosInPartitionTable);

            oldPartitionTableSize = newPartitionTableSize;
            newPartitionTableSize = this->lowValTable.size();
            lastPosInPartitionTable += newPartitionTableSize - oldPartitionTableSize;
            pair = this->crackForHighValue(tempLowVal, lastPosInPartitionTable, highVal, false);
            long highValPosInArray = pair.second;

            return {lowValPosInArray, highValPosInArray};
        }


        
        std::pair<rs::Coord<KeyType>*, size_t>  buildIndexForGap(long firstPosInPartitionTable, long lastPosInPartitionTable)
        {
            rs::Coord<KeyType>* res = nullptr;
            size_t resSize = 0;

            for (long i = firstPosInPartitionTable ; i < lastPosInPartitionTable ; i++)
            {
                long startPos = this->highPosTable[i] + 1, endPos = this->lowPosTable[i+1];
                if(startPos != endPos)
                {
                    std::pair<KeyType, KeyType> p = this->getMinMaxElement(startPos, endPos-1);
                    KeyType lowVal = p.first, highVal = p.second;
                    std::tuple<long, long, rs::Coord<KeyType>*, size_t> p = this->buildLearnedIndex(lowVal, highVal, startPos, endPos-1);
                    lastPosInPartitionTable += 1;

                    rs::Coord<KeyType>* insertedData = std::get<2>(p);
                    size_t size = std::get<3>(p);

                    rs::Coord<KeyType>* temp = new rs::Coord<KeyType> [resSize + size];
                    std::copy(res, res+resSize, temp);
                    std::copy(insertedData, insertedData+size,temp+resSize);
                    resSize += size;
                    res = temp;
                    //this->queryTypeCount[BUILD_LEARNED_INDEX]++; 
                }
            }
            return {res, resSize};
            
        }
#else

        std::tuple<long, long, rs::Coord<KeyType>*, size_t>  executeOverlapQuery(KeyType lowVal, KeyType highVal)
        {
            long firstPosInPartitionTable = this->getFirstEntryGreaterThanLowVal(lowVal);
            long lastPosInPartitionTable = this->getLastEntryLessThanHighVal(highVal);
            
            std::tuple<long, long, rs::Coord<KeyType>*, size_t> firstResults = crackForLowValue(lowVal, highValTable[lastPosInPartitionTable], lastPosInPartitionTable);
            long lowValPosInArray = std::get<0>(firstResults);
            rs::Coord<KeyType>* firstInserted = std::get<2>(firstResults);
            size_t firstSize = std::get<3>(firstResults);

            long startPos = this->highPosTable[lastPosInPartitionTable] + 1;
            long endPos = (lastPosInPartitionTable < this->lowPosTable.size() - 1) ? 
                            this->lowPosTable[lastPosInPartitionTable + 1] - 1 : this->data.size() - 2;

            KeyType startVal = this->getMinElement(startPos, endPos);
            //KeyType startVal = *std::min_element(data.begin() + startPos, data.begin() + endPos + 1);
            std::tuple<long, long, rs::Coord<KeyType>*, size_t> secondResults = this->buildLearnedIndex(startVal, highVal, startPos, endPos);
            long highValPosInArray = std::get<1>(secondResults);
            rs::Coord<KeyType>* secondInserted = std::get<2>(secondResults);
            size_t secondSize = std::get<3>(secondResults);

            size_t finalInsertedAndPendingResultsSize = firstSize + secondSize;
            rs::Coord<KeyType>* finalInsertedAndPendingResults = new rs::Coord<KeyType> [finalInsertedAndPendingResultsSize];
            std::copy(firstInserted, firstInserted+firstSize, finalInsertedAndPendingResults);
            std::copy(secondInserted, secondInserted+secondSize, finalInsertedAndPendingResults+firstSize);

            return {lowValPosInArray, highValPosInArray, finalInsertedAndPendingResults, finalInsertedAndPendingResultsSize};
        }
#endif


        
        std::pair<rs::Coord<KeyType>*, size_t>  buildIndexForGap(long firstPosInPartitionTable, long lastPosInPartitionTable)
        {
            rs::Coord<KeyType>* res = nullptr;
            size_t resSize = 0;

            for (long i = firstPosInPartitionTable ; i < lastPosInPartitionTable ; i++)
            {
                long startPos = this->highPosTable[i] + 1, endPos = this->lowPosTable[i+1];
                if(startPos != endPos)
                {
                    std::pair<KeyType, KeyType> p = this->getMinMaxElement(startPos, endPos-1);
                    KeyType lowVal = p.first, highVal = p.second;
                    std::tuple<long, long, rs::Coord<KeyType>*, size_t> tup = this->buildLearnedIndex(lowVal, highVal, startPos, endPos-1);
                    lastPosInPartitionTable += 1;

                    rs::Coord<KeyType>* insertedData = std::get<2>(tup);
                    size_t size = std::get<3>(tup);

                    rs::Coord<KeyType>* temp = new rs::Coord<KeyType> [resSize + size];
                    std::copy(res, res+resSize, temp);
                    std::copy(insertedData, insertedData+size,temp+resSize);
                    resSize += size;
                    res = temp;
                    //this->queryTypeCount[BUILD_LEARNED_INDEX]++; 
                }
            }
            return {res, resSize};
            
        }

        bool checkAnswers(KeyType lowVal, KeyType highVal, long lowValPos, long highValPos)
        {
            auto startIter = this->data.begin() + lowValPos, endIter = this->data.begin() + highValPos + 1;
            return (this->data[lowValPos].x == lowVal && this->data[highValPos].x == highVal &&
                std::is_sorted(startIter, endIter, compare_by_key));
        }

        void getQueryStats()
        {
            for (int i = 0 ; i < MAX_QUERY_TYPE ; i++)
            {
                std::cout << "Query Type " << i << " = " << this->queryTypeCount[i] << " Runtime = " << this->queryRunTime[i] << std::endl;            
            }
            //std::cout << "Partition Time = " << this->partitionTime << endl;
            //std::cout << "Sort Time = " << this->sortTime << endl;
            //std::cout << "Learned Index Build Time = " << this->learnedIndexBuildTime << endl;
        }


};

#endif
