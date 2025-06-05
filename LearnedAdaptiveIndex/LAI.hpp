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

#define SORTING_THRESHOLD 6000


template <typename KeyType>
class LAI{

    private:
        KeyType* data;
        int dataset_size;
        std::vector<KeyType> lowValTable;
        std::vector<KeyType> highValTable;
        std::vector<long> lowPosTable;
        std::vector<long> highPosTable;
        std::vector<baseLearnedModel<KeyType>* > learnedIndexTable;
        int error;
        queryType qType;
        bool useLearnedSort;
        unsigned long queryTypeCount[MAX_QUERY_TYPE];
        double queryRunTime[MAX_QUERY_TYPE];
        //double partitionTime, sortTime, learnedIndexBuildTime;
          
    public:

        LAI(KeyType* data, int size, int error, bool useLearnedSort=true)
        {
            std::cout << "LAI with array" << std::endl;
            this->data = data;
            this->dataset_size = size;
            this->data[0] = -INFINITY;
            this->data[size+1] = INFINITY;
            this->error = error;
            this->useLearnedSort = useLearnedSort;
            //this->partitionTime = this->sortTime = this->learnedIndexBuildTime = 0;

            for (int i = 0 ; i < MAX_QUERY_TYPE ; i++)
            {
                this->queryTypeCount[i] = 0;
                this->queryRunTime[i] = 0;
            }
            int reserve_size = 20000;
            this->lowValTable.reserve(reserve_size);
            this->highValTable.reserve(reserve_size);
            this->lowPosTable.reserve(reserve_size);
            this->highPosTable.reserve(reserve_size);
            this->learnedIndexTable.reserve(reserve_size);
        }


        std::pair<long, long> query(KeyType lowVal, KeyType highVal, bool measure_stats=true)
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
            std::pair<long, long> pos;
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
            

            return pos;
        }


 #if 0       
        long partition(long startPos, long endPos, KeyType pivot)
        {    
            long long pivotPos;

            //find whether the pivot is present in the sub-array or not
            //auto startIter(data.begin() + startPos);
            //auto endIter(data.begin() + endPos+1);
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
                    if (data[i] > pivot)
                        break;
                    if (data[i] == pivot)
                        pivotPos = i;    
                }                        
                while(true)
                {
                    j-= 1;
                    if (data[j] <= pivot)
                        break;
                }        

                if (data[j] == pivot)
                    pivotPos = j;

                if (i < j)
                { 
                    std::swap(data[i], data[j]);
                       
                    if (data[i] == pivot)
                        pivotPos = i;
                    else if (data[j] == pivot)
                        pivotPos = j; 
                       
                }   

            }
            if (found)
            {
                std::swap(data[j], data[pivotPos]);
            }
                
            return j;
        }    
#endif

#if 0
    long partition(long startPos, long endPos, KeyType pivot)
        {    
            
            auto middle1 = std::partition(data+startPos, data+endPos+1, [pivot](KeyType& em)
                {
                    return em < pivot;
                });
            auto middle2 = std::partition(middle1, data+endPos+1, [pivot](KeyType& em)
                {
                    //return !(pivot < em);
                    return pivot == em;
                });
            
            long pivotPos = middle2 - data - 1;
            return  pivotPos;              
        }    
#endif



        std::pair<long, long>  crackAndSort(KeyType lowVal, KeyType highVal, long startPos, long endPos)
        {
            if (startPos == -1 and endPos == -1) 
            {
                std::pair<long, long> p = getPositions(lowVal, highVal);
                startPos = p.first;
                endPos = p.second;
            }

            //auto started = std::chrono::high_resolution_clock::now();
            //long p = this->partition(startPos, endPos, lowVal);
            //long q = this->partition(p+1, endPos, highVal);

            #if 1
            auto m1 = std::partition(data+startPos, data+endPos+1, [lowVal](KeyType& em)
                {
                    return em < lowVal;
                });
            auto m2 = std::partition(m1, data+endPos+1, [highVal](KeyType& em)
                {
                    return em < highVal;
                });
            auto m3 = std::partition(m2, data+endPos+1, [highVal](KeyType& em)
                {
                    return em == highVal;
                });
            auto m4 = std::partition(m1, m3, [lowVal](KeyType& em)
                {
                    return em == lowVal;
                });

            long p = m4 - data - 1, q = m3 - data - 1;    
            #endif        



            if (data[p] != lowVal)
                p += 1;
        
            if (p < q)
            {
                auto startIter = data + p;
                auto endIter = data + q + 1;
                long numElementsToSort = q - p + 1;
                //started = std::chrono::high_resolution_clock::now(); 
                
                #if 0
                if(!this->useLearnedSort)
                    std::sort(startIter, endIter);
                else
                    learnedSort(p, q);
                #endif    
                    
                #if 1
                if((numElementsToSort > SORTING_THRESHOLD) /*&& this->useLearnedSort*/)
                    learnedSort(p, q);
                else
                   std::sort(startIter, endIter);
                #endif        
            }
            return std::pair<long, long>(p,q);
        }
                
            
        void learnedSort(long startPos, long endPos)
        {
            #ifdef DEBUG_LEARNED_SORT
            vector<KeyType> data_2(data, data + size + 2);
            auto tempIter1 = data_2.begin() + startPos + 1, tempIter2 = data_2.begin() + endPos;
            auto started = std::chrono::high_resolution_clock::now(); 
            std::sort(tempIter1, tempIter2);
            auto done = std::chrono::high_resolution_clock::now();
            double timeToNativeSort = (double)std::chrono::duration_cast<std::chrono::nanoseconds>(done-started).count();
            
            started = std::chrono::high_resolution_clock::now();
            #endif 
            
            KeyType startVal = data[startPos], endVal = data[endPos];
            double slope = (double)(endPos - startPos) / (double)(endVal - startVal);
            double intercept = endPos - (slope * endVal);

            vector<KeyType> overflow_buffer;
            long size = endPos - startPos - 1;
            overflow_buffer.reserve(size);
            overflow_buffer.clear();

            //KeyType maxVal = data[data.size() - 1];
            KeyType maxVal = data[dataset_size + 1];
            vector<KeyType> temp(size, maxVal);

            //The elements at startPos and endPos are already in their correct locations, so no need
            //to consider them
            for (long i = startPos + 1 ; i < endPos ; i++)
            {
                KeyType key = data[i];
                long pos = round(slope * key + intercept) - (startPos + 1);
                if (temp[pos] == maxVal)
                    temp[pos] = key;
                else
                {
                    overflow_buffer.push_back(key);
                }
            }

            temp.erase(std::remove(temp.begin(), temp.end(), maxVal), temp.end());
            std::sort(overflow_buffer.begin(), overflow_buffer.end());
            vector<KeyType> temp2(size);
            
            KeyType* startIter = data + startPos + 1;

            std::merge(temp.begin(), temp.end(), overflow_buffer.begin(), overflow_buffer.end(), startIter);
            
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
                //pos = this->data.size() - 2;
                pos = dataset_size;
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
                endPos = dataset_size;
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
            KeyType max_elem = this->data[start];
            for (int i = start+1 ; i<= end ; i++)
            {
                KeyType elem = this->data[i];
                max_elem = std::max(max_elem, this->data[i]);
            }
            return max_elem;
        }


        
        KeyType  getMinElement(long start, long end)
        {
            KeyType min_elem = this->data[start];
            for (int i = start+1 ; i<= end ; i++)
            {
                KeyType elem = this->data[i];
                min_elem = std::min(min_elem, this->data[i]);
            }
            return min_elem;
        }


        
        std::pair<KeyType, KeyType>  getMinMaxElement(long start, long end)
        {
            KeyType min_elem = this->data[start], max_elem = this->data[start];
            for (int i = start+1 ; i<= end ; i++)
            {
                KeyType elem = this->data[i];
                min_elem = std::min(min_elem, this->data[i]);
                max_elem = std::max(max_elem, this->data[i]);
            }
            return {min_elem, max_elem};
        }


        
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
            learnedIndexModel->buildIndex(data, startPos, endPos);
            return learnedIndexModel;
        } 


        
        std::pair<long, long>  buildLearnedIndex(KeyType lowVal, KeyType highVal, long lowValPos=-1, long highValPos=-1)
        {
            std::pair<long, long> pos = crackAndSort(lowVal, highVal, lowValPos, highValPos);
            long startPos = pos.first, endPos = pos.second;
            bool is_sorted = std::is_sorted(this->data+startPos, this->data+endPos+1);
            
            if (startPos > endPos)
            {
                cout << "No element in this range" << endl;
                return {-1, -1};
            }
                
            if ((startPos == endPos) && 
                    ((this->qType == CRACK_LOW_VALUE) || (this->qType == CRACK_HIGH_VALUE)
                    || (this->qType == BUILD_LEARNED_INDEX) || (this->qType == OVERLAP)))
                return {startPos, endPos};
                
            baseLearnedModel<KeyType> *learnedIndexModel = this->buildLearnedIndexModel(startPos, endPos);

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


        
        long  getFirstEntryGreaterThanLowVal(KeyType lowVal)
        {
            auto insertPosIter = std::upper_bound (this->lowValTable.begin(), this->lowValTable.end(), lowVal);
            long pos = insertPosIter -  this->lowValTable.begin();
            return pos;                                         
        }


        
        long  getLastEntryLessThanHighVal(KeyType highVal)
        {
            auto insertPosIter = std::upper_bound (this->lowValTable.begin(), this->lowValTable.end(), highVal);
            long pos = insertPosIter -  this->lowValTable.begin() - 1;
            return pos;                                        
        }


        
        std::pair<long, long>  getResultsFromSameBound(KeyType lowVal, KeyType highVal, long posInPartitionTable)
        {
            baseLearnedModel<KeyType> *learnedIndexModel = this->learnedIndexTable[posInPartitionTable];
            long lowValPosInArray = learnedIndexModel->getPosition(lowVal);
            long highValPosInArray = learnedIndexModel->getPosition(highVal);
            return {lowValPosInArray, highValPosInArray};
        }


        
        std::pair<long, long>  getResultsFromDifferentBounds(KeyType lowVal, KeyType highVal, long lowValPosInPartitionTable, long highValPosInPartitionTable)
        {
            baseLearnedModel<KeyType> *learnedIndexModel = this->learnedIndexTable[lowValPosInPartitionTable];
            long lowValPosInArray = learnedIndexModel->getPosition(lowVal);
            
            learnedIndexModel = this->learnedIndexTable[highValPosInPartitionTable];
            long highValPosInArray = learnedIndexModel->getPosition(highVal);

            this->buildIndexForGap(lowValPosInPartitionTable, highValPosInPartitionTable);
            return {lowValPosInArray, highValPosInArray};
        }


        
        std::pair<long, long>  crackForLowValue(KeyType lowVal, KeyType highVal, long highValPosInPartitionTable, bool computeHigh=true)
        {
            long highValPosInArray = -1;
            if (computeHigh)
            {
                baseLearnedModel<KeyType> *learnedIndexModel = this->learnedIndexTable[highValPosInPartitionTable];
                highValPosInArray = learnedIndexModel->getPosition(highVal);
            }

            long firstPosInPartitionTable = this->getFirstEntryGreaterThanLowVal(lowVal);
            this->buildIndexForGap(firstPosInPartitionTable, highValPosInPartitionTable);

            long startPos = (firstPosInPartitionTable == 0) ? 1 : this->highPosTable[firstPosInPartitionTable - 1] + 1;
            long endPos = this->lowPosTable[firstPosInPartitionTable] - 1;

            KeyType endVal = this->getMaxElement(startPos, endPos);
            std::pair p = this->buildLearnedIndex(lowVal, endVal, startPos, endPos);
            long lowValPosInArray = p.first;

            //if one point is left to be cracked, add that point in the beginning and rebuild the whole learned model again
            if (lowVal == endVal)
            {
                
                baseLearnedModel<KeyType> *learnedIndexModel = this->buildLearnedIndexModel(lowValPosInArray, this->highPosTable[firstPosInPartitionTable]);
                this->lowValTable[firstPosInPartitionTable] = lowVal;
                this->lowPosTable[firstPosInPartitionTable] = lowValPosInArray;
                this->learnedIndexTable[firstPosInPartitionTable] = learnedIndexModel;
            }

            return {lowValPosInArray, highValPosInArray};
        }


        
        std::pair<long, long>  crackForHighValue(KeyType lowVal, long lowValPosInPartitionTable, KeyType highVal, bool computeLow=true)
        {
            long lowValPosInArray = -1;
            if (computeLow)
            {
                baseLearnedModel<KeyType> *learnedIndexModel = this->learnedIndexTable[lowValPosInPartitionTable];
                lowValPosInArray = learnedIndexModel->getPosition(lowVal);
            }

            long lastPosInPartitionTable = this->getLastEntryLessThanHighVal(highVal);
            long oldPartitionTableSize = this->lowValTable.size();
            this->buildIndexForGap(lowValPosInPartitionTable, lastPosInPartitionTable);
            long newPartitionTableSize = this->lowValTable.size();
            lastPosInPartitionTable += newPartitionTableSize - oldPartitionTableSize; 

            long startPos = this->highPosTable[lastPosInPartitionTable] + 1;
            long endPos = (lastPosInPartitionTable < this->lowPosTable.size() - 1) ? 
                            this->lowPosTable[lastPosInPartitionTable + 1] - 1 : this->dataset_size;

            KeyType startVal = this->getMinElement(startPos, endPos);
            std:: pair p = this->buildLearnedIndex(startVal, highVal, startPos, endPos);
            long highValPosInArray = p.second;
            
            //if one point is left to be cracked, add that point to the existing partition info
            if (startVal == highVal)
            {
                baseLearnedModel<KeyType> *learnedIndexModel = this->buildLearnedIndexModel(this->lowPosTable[lastPosInPartitionTable], highValPosInArray);
                this->highValTable[lastPosInPartitionTable] = highVal;
                this->highPosTable[lastPosInPartitionTable] = highValPosInArray;
                this->learnedIndexTable[lastPosInPartitionTable] = learnedIndexModel;
            }
                           
            return {lowValPosInArray, highValPosInArray};
        }


        
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


        
        void  buildIndexForGap(long firstPosInPartitionTable, long lastPosInPartitionTable)
        {
            for (long i = firstPosInPartitionTable ; i < lastPosInPartitionTable ; i++)
            {
                long startPos = this->highPosTable[i] + 1, endPos = this->lowPosTable[i+1];
                if(startPos != endPos)
                {
                    std::pair<KeyType, KeyType> p = this->getMinMaxElement(startPos, endPos-1);
                    KeyType lowVal = p.first, highVal = p.second;
                    this->buildLearnedIndex(lowVal, highVal, startPos, endPos-1);
                    lastPosInPartitionTable += 1;
                }
            }
            
        }


        bool checkAnswers(KeyType lowVal, KeyType highVal, long lowValPos, long highValPos)
        {
            KeyType *startIter = this->data + lowValPos, *endIter = this->data + highValPos + 1;
            return (this->data[lowValPos] == lowVal && this->data[highValPos] == highVal &&
                std::is_sorted(startIter, endIter));
        }

        void getQueryStats()
        {
            for (int i = 0 ; i < MAX_QUERY_TYPE ; i++)
            {
                std::cout << "Query Type " << i << " = " << this->queryTypeCount[i] << " Runtime = " << this->queryRunTime[i] << std::endl;            
            }
        }


};

#endif
