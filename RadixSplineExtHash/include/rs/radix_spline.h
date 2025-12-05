#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>

#include "common.h"
#include "ext_hash.h"

namespace rs {

// Approximates a cumulative distribution function (CDF) using spline
// interpolation.
template <class KeyType>
class RadixSpline {
 public:
  RadixSpline() = default;

  RadixSpline(KeyType min_key, KeyType max_key, size_t num_keys,
              size_t num_radix_bits, size_t num_shift_bits, size_t max_error,
              std::vector<uint32_t> radix_table,
              rs::Directory<KeyType>spline_points)
      : min_key_(min_key),
        max_key_(max_key),
        num_keys_(num_keys),
        num_radix_bits_(num_radix_bits),
        num_shift_bits_(num_shift_bits),
        max_error_(max_error),
        radix_table_(std::move(radix_table)),
        spline_points_(std::move(spline_points)),
        overflow_buffer_size(spline_points.getOverflowBufferSize())
        {
        }

  // Returns the estimated position of `key`.
  double GetEstimatedPosition(const Coord<KeyType> down, const Coord<KeyType> up, const KeyType key) const {
    // Truncate to data boundaries.
    if (key <= min_key_) return 0;
    if (key >= max_key_) return num_keys_ - 1;

    
    // Find spline segment with `key` ∈ (spline[index - 1], spline[index]].
    /*  //skd
    const size_t index = GetSplineSegment(key);
    const Coord<KeyType> down = spline_points_[index - 1];
    const Coord<KeyType> up = spline_points_[index];
    */
    //const Coord<KeyType> down = p.first;
    //const Coord<KeyType> up = p.second;

    // Compute slope.
    const double x_diff = up.x - down.x;
    const double y_diff = up.y - down.y;
    const double slope = y_diff / x_diff;

    // Interpolate.
    const double key_diff = key - down.x;
    return std::fma(key_diff, slope, down.y);
  }

  // Returns a search bound [begin, end) around the estimated position.
  SearchBound GetSearchBound(const size_t estimate) const {
    //const size_t estimate = GetEstimatedPosition(key);
    const size_t begin = (estimate < max_error_) ? 0 : (estimate - max_error_);
    // `end` is exclusive.
    const size_t end = (estimate + max_error_ + 2 > num_keys_)
                           ? num_keys_
                           : (estimate + max_error_ + 2);
    return SearchBound{begin, end};
  }

  size_t getPositionGreaterThanOrEquals(std::vector<Coord<KeyType>> &keys, size_t startPos, const KeyType search_key)
  {
    std::pair<rs::BucketEntry<KeyType>*, rs::BucketEntry<KeyType>*> p = GetSplineSegment(search_key);
    rs::BucketEntry<KeyType>* first = p.first;
    rs::BucketEntry<KeyType>* second = p.second;
    Coord<KeyType> c1 = {first->val, first->pos}, c2 = {second->val, second->pos};
    size_t estimatedPosition = GetEstimatedPosition(c1, c2, search_key);
    SearchBound bound = GetSearchBound(estimatedPosition);
    
    //Thread1
    auto start = begin(keys) + startPos + bound.begin, last = begin(keys) + startPos + bound.end;
    Coord<KeyType> c = {search_key, -1};
    size_t posInArray = std::lower_bound(start, last, c, compare_by_key) - begin(keys);
    //if (keys[posInArray].x == search_key)
      return posInArray;
    //else //(keys[posInArray].x > search_key)
      //return --posInArray;
    
  }


  size_t getPosition(std::vector<Coord<KeyType>> &keys, size_t startPos, const KeyType search_key)
  {
    std::pair<rs::BucketEntry<KeyType>*, rs::BucketEntry<KeyType>*> p = GetSplineSegment(search_key);
    rs::BucketEntry<KeyType>* first = p.first;
    rs::BucketEntry<KeyType>* second = p.second;
    Coord<KeyType> c1 = {first->val, first->pos}, c2 = {second->val, second->pos};
    size_t estimatedPosition = GetEstimatedPosition(c1, c2, search_key);
    SearchBound bound = GetSearchBound(estimatedPosition);
    
    //Thread1
    auto start = begin(keys) + startPos + bound.begin, last = begin(keys) + startPos + bound.end;
    Coord<KeyType> c = {search_key, -1};
    size_t posInArray = std::lower_bound(start, last, c, compare_by_key) - begin(keys);
    if (keys[posInArray].x == search_key)
      return posInArray;
    //else //(keys[posInArray].x > search_key)
      //return --posInArray;

    
    //Thread2
    Coord<KeyType> *overflow_buffer;
    size_t size;
    if (spline_points_.getPrefix(first->val) == spline_points_.getPrefix(second->val))
    {
      overflow_buffer = second->left_overflow_buffer;
      size = second->getLeftBufferSize();
    }
    else
    {
      Coord<KeyType>*temp_buffer1 = first->right_overflow_buffer;
      size_t size1 = first->getRightBufferSize();
      Coord<KeyType>*temp_buffer2 = second->left_overflow_buffer;
      size_t size2 = second->getLeftBufferSize();
      size = size1 + size2;
      
      overflow_buffer = new Coord<KeyType> [size];
      std::copy(temp_buffer1, temp_buffer1 + size1, overflow_buffer);
      std::copy(temp_buffer2, temp_buffer2 + size2, overflow_buffer + size1);
    }

    double posInBuffer = -1;
    for(int i = 0; i < size; i++)
      if (overflow_buffer[i].x == search_key)
      {
        posInBuffer = overflow_buffer[i].y;
        break;
      }

    return posInBuffer;
  }

  std::pair<Coord<KeyType>*, size_t> RangeSearch(std::vector<Coord<KeyType>> &keys, const KeyType low, const KeyType high, bool search_inserted_data_only=false)
  {
    /* Start searching in the original data array followed by the overflow buffer*/
    
    //Search for low in the index
    size_t posInArrayLow, posInArrayHigh;
    rs::BucketEntry<KeyType> *first_low, *second_low, *first_high, *second_high;
    if (low < min_key_)
    {
      first_low = second_low = spline_points_.getBucket(radix_table_[0])->getElement(0);
      posInArrayLow = first_low->pos;
    }
    else
    {
      std::pair<rs::BucketEntry<KeyType>*, rs::BucketEntry<KeyType>*> p = GetSplineSegment(low);
      first_low = p.first;
      second_low = p.second;
      if (!search_inserted_data_only)
      {
        Coord<KeyType> c1 = {first_low->val, first_low->pos}, c2 = {second_low->val, second_low->pos};
        size_t estimatedPosition = GetEstimatedPosition(c1, c2, low);
        SearchBound bound = GetSearchBound(estimatedPosition);
        
        auto start = begin(keys) + bound.begin, last = begin(keys) + bound.end;
        Coord<KeyType> c = {low, -1};
        posInArrayLow = std::lower_bound(start, last, c, compare_by_key) - begin(keys);
      }
    }

    //Search for high in the index
    if (high > max_key_)
    {
      first_high = second_high = spline_points_.getBucket(radix_table_.back())->getLastElement();
      posInArrayHigh = first_high->pos;
    }
    else
    {
      std::pair<rs::BucketEntry<KeyType>*, rs::BucketEntry<KeyType>*> p = GetSplineSegment(high);
      first_high = p.first;
      second_high = p.second;
      if (!search_inserted_data_only)
      {
        Coord<KeyType> c1 = {first_high->val, first_high->pos}, c2 = {second_high->val, second_high->pos};
        size_t estimatedPosition = GetEstimatedPosition(c1, c2, high);
        SearchBound bound = GetSearchBound(estimatedPosition);
        
        auto start = begin(keys) + bound.begin, last = begin(keys) + bound.end;
        Coord<KeyType> c = {high, -1};
        posInArrayHigh = std::lower_bound(start, last, c, compare_by_key) - begin(keys);
        if (keys[posInArrayHigh].x > high)
          posInArrayHigh--;
      }
    }

    size_t resultSizeInArray = 0;
    Coord<KeyType> *resInArray = nullptr;
    if (!search_inserted_data_only)
    {
      resultSizeInArray = posInArrayHigh - posInArrayLow + 1;
      resInArray = new Coord<KeyType> [resultSizeInArray];
      std::copy(keys.begin() + posInArrayLow, keys.begin() + posInArrayHigh + 1, resInArray); 
    } 

  
    //Search in the overflow buffers associated with all the bucket entries in-between
    KeyType first_low_prefix = spline_points_.getPrefix(first_low->val);
    KeyType second_low_prefix = spline_points_.getPrefix(second_low->val);
    KeyType first_high_prefix = spline_points_.getPrefix(first_high->val);
    KeyType second_high_prefix = spline_points_.getPrefix(second_high->val);
    auto it = std::lower_bound(radix_table_.begin(), radix_table_.end(), first_low_prefix);
    size_t posInRadixTable = it - radix_table_.begin();
    std::vector<Coord<KeyType>> resInBuffer;
    resInBuffer.reserve(keys.size());

    rs::Coord<KeyType> *buffer;
    size_t size;

    if (false && (first_low_prefix == second_high_prefix)) // if both low and high belong to same bucket
    {
       //std::cout << "This case is not handled yet" << std::endl;
       //std::cout << __FILE__ << " " << __LINE__ << std::endl; 
        buffer = first_low->left_overflow_buffer;
        size = first_low->getLeftBufferSize();
        for (size_t i = 0 ; i < size ; i++)
        {
          if (buffer[i].x >= low && buffer[i].x <= high)
            resInBuffer.push_back(buffer[i]);
        }

        buffer = first_low->right_overflow_buffer;
        size = first_low->getRightBufferSize();
        for (size_t i = 0 ; i < size ; i++)
        {
          if (buffer[i].x >= low && buffer[i].x <= high)
            resInBuffer.push_back(buffer[i]);
        }
    }

    else
    {



      rs::BucketEntry<KeyType> *currentBucketEntry;
      if (first_low_prefix != second_low_prefix)
      {
        buffer = first_low->right_overflow_buffer;
        size = first_low->getRightBufferSize();
        for (size_t i = 0 ; i < size ; i++)
        {
          if (buffer[i].x >= low && buffer[i].x <= high)
            resInBuffer.push_back(buffer[i]);
        }
        rs::Bucket<KeyType> *current_bucket = spline_points_.getBucket(radix_table_[posInRadixTable]);
        rs::Bucket<KeyType> *next_bucket = spline_points_.getBucket(radix_table_[++posInRadixTable]);
        while (current_bucket == next_bucket)
        {
          posInRadixTable++;
          current_bucket = next_bucket;
          next_bucket = spline_points_.getBucket(radix_table_[posInRadixTable+1]);
        }
      }
      buffer = second_low->left_overflow_buffer;
      size = second_low->getLeftBufferSize();
      for (size_t i = 0 ; i < size ; i++)
      {
        if (buffer[i].x >= low && buffer[i].x <= high)
          resInBuffer.push_back(buffer[i]);
      }

      buffer = second_low->right_overflow_buffer;
      size = second_low->getRightBufferSize();
      if(second_low->val < first_high->val)
        resInBuffer.insert(resInBuffer.end(), buffer, buffer+size);    
      else if (second_low->val == first_high->val)
      {
        for (size_t i = 0 ; i < size ; i++)
        {
          if (buffer[i].x >= low && buffer[i].x <= high)
            resInBuffer.push_back(buffer[i]);
        }
      } 
      currentBucketEntry = second_low + 1;

      while(radix_table_[posInRadixTable] < second_high_prefix)
      {
        KeyType current_prefix = radix_table_[posInRadixTable];
        rs::BucketEntry<KeyType> *lastBucketEntry = spline_points_.getBucket(current_prefix)->getLastElement();
        while(currentBucketEntry <= lastBucketEntry)
        {
          buffer = currentBucketEntry->left_overflow_buffer;
          size = currentBucketEntry->getLeftBufferSize();
          resInBuffer.insert(resInBuffer.end(), buffer, buffer+size);
  
          buffer = currentBucketEntry->right_overflow_buffer;
          size = currentBucketEntry->getRightBufferSize();
          if (currentBucketEntry == first_high)
          {
            for (size_t i = 0 ; i < size ; i++)
            {
              if (buffer[i].x <= high)
                resInBuffer.push_back(buffer[i]);
            }
          }
          else
            resInBuffer.insert(resInBuffer.end(), buffer, buffer+size);

          currentBucketEntry++;
        }

        if (radix_table_[posInRadixTable + 1] == second_high_prefix)
        {
          posInRadixTable++;
          currentBucketEntry = spline_points_.getBucket(radix_table_[posInRadixTable])->getElement(0);
          break;
        }

        rs::Bucket<KeyType> *current_bucket = spline_points_.getBucket(radix_table_[posInRadixTable]);
        rs::Bucket<KeyType> *next_bucket = spline_points_.getBucket(radix_table_[++posInRadixTable]);
        while (current_bucket == next_bucket)
        {
          posInRadixTable++;
          current_bucket = next_bucket;
          next_bucket = spline_points_.getBucket(radix_table_[posInRadixTable]);
        }

        currentBucketEntry = spline_points_.getBucket(radix_table_[posInRadixTable])->getElement(0);
      }
      
      //When code reaches here, currentBucketEntry points to the first entry of the end_prefix.
      //Copy all left_buffers until currentBucketEntry == second_high
      while(currentBucketEntry < second_high)
      {
          buffer = currentBucketEntry->left_overflow_buffer;
          size = currentBucketEntry->getLeftBufferSize();
          resInBuffer.insert(resInBuffer.end(), buffer, buffer+size);;
          currentBucketEntry++;
      }

      //buffer = currentBucketEntry->left_overflow_buffer;
      if (second_low != second_high)
      {
        buffer = second_high->left_overflow_buffer;
        size = second_high->getLeftBufferSize();
        for (size_t i = 0 ; i < size ; i++)
        {
          if (buffer[i].x <= high)
            resInBuffer.push_back(buffer[i]);
        }
      }

      if (high > max_key_)
      {
        buffer = second_high->right_overflow_buffer;
        size = second_high->getRightBufferSize();
        for (size_t i = 0 ; i < size ; i++)
        {
          if (buffer[i].x <= high)
            resInBuffer.push_back(buffer[i]);
        }      
      }
    }
    
    if (!search_inserted_data_only)
    {
      size_t total_result_size = resultSizeInArray + resInBuffer.size();
      Coord<KeyType> *final_result = new Coord<KeyType> [total_result_size];
      std::copy(resInArray, resInArray+resultSizeInArray, final_result);
      std::copy(resInBuffer.begin(), resInBuffer.end(), final_result+resultSizeInArray);
      return std::pair<Coord<KeyType>*, size_t> (final_result, total_result_size);
    }
    else
    {
      size_t total_result_size = resInBuffer.size();
      Coord<KeyType> *final_result = new Coord<KeyType> [total_result_size];
      std::copy(resInBuffer.begin(), resInBuffer.end(), final_result);
      return std::pair<Coord<KeyType>*, size_t> (final_result, total_result_size);
    }

  }

  // Returns the size in bytes.
  /*  //skd
  size_t GetSize() const {
    return sizeof(*this) + radix_table_.size() * sizeof(uint32_t) +
           spline_points_.size() * sizeof(Coord<KeyType>);
  }
  */

  size_t GetSize() const {
    return sizeof(*this) + radix_table_.size() * sizeof(uint32_t) +
           spline_points_.globalElementCount() * sizeof(Coord<KeyType>);
  }

  std::pair<Coord<KeyType>*, size_t> getAllElementsLessThanMin()
  {
    KeyType prefix = radix_table_[0];
    rs::Bucket<KeyType>* bucket = spline_points_.getBucket(prefix);
    rs::BucketEntry<KeyType>* entry = bucket->getElement(0);
    Coord<KeyType>* res = entry->left_overflow_buffer;
    size_t size = entry->getLeftBufferSize();
    return {res, size};
  }


  std::pair<Coord<KeyType>*, size_t> getAllElementsGreaterThanMax()
  {
    KeyType prefix = radix_table_.back();
    rs::Bucket<KeyType>* bucket = spline_points_.getBucket(prefix);
    rs::BucketEntry<KeyType>* entry = bucket->getLastElement();
    Coord<KeyType>* res = entry->right_overflow_buffer;
    size_t size = entry->getRightBufferSize();
    return {res, size};
  }

  int insertNewKey(KeyType newKey, double pos)
  {
    KeyType prefix;
    rs::Bucket<KeyType>* bucket;
    
    if (newKey < min_key_)
    {
      prefix = radix_table_[0];
      bucket = spline_points_.getBucket(prefix);
      rs::BucketEntry<KeyType>* entry = bucket->getElement(0);
      if (entry->getLeftBufferSize() == overflow_buffer_size)
        return retrain();
      else
      {
        entry->insertLeft(newKey, pos); 
        return 1;
      } 
    }
    
    if (newKey > max_key_)
    {
      prefix = radix_table_.back();
      bucket = spline_points_.getBucket(prefix);
      rs::BucketEntry<KeyType>* entry = bucket->getLastElement();
      if (entry->getRightBufferSize() == overflow_buffer_size)
        return retrain();
      else
      {
        entry->insertRight(newKey, pos); 
        return 1;
      } 
    }

    prefix = spline_points_.getPrefix(newKey);
    bucket = spline_points_.getBucket(prefix);
    
    if (bucket == nullptr || bucket->getNumElements() == 0) //if new key belongs to an empty bucket
    {
      auto it = std::lower_bound(radix_table_.begin(), radix_table_.end(), prefix);
      if (it > radix_table_.begin() && it < radix_table_.end()) //new key falls within this bucket
      {
        bucket = spline_points_.getBucket(*it);
        rs::BucketEntry<KeyType>* right = bucket->getElement(0); //first entry of the right bucket

        bucket = spline_points_.getBucket(*(--it));
        rs::BucketEntry<KeyType>* left = bucket->getLastElement(); //last entry of left bucket

        size_t left_entry_right_buffer_size = left->getRightBufferSize(), right_entry_left_buffer_size = right->getLeftBufferSize();
        if (left_entry_right_buffer_size == overflow_buffer_size && right_entry_left_buffer_size == overflow_buffer_size)
        {
          return retrain();
        }
        
        if (left_entry_right_buffer_size <= right_entry_left_buffer_size)
        {
          left->insertRight(newKey, pos);
        }
        else
        {
          right->insertLeft(newKey, pos);
        }
      }
      else if(it == radix_table_.begin()) //new key falls before the first bucket
      {
          bucket = spline_points_.getBucket(*it);
          rs::BucketEntry<KeyType>* right = bucket->getElement(0); //first entry of the right bucket
          if (right->getLeftBufferSize() == overflow_buffer_size)
            return retrain();
          else
            right->insertLeft(newKey, pos);
      }
      else  //new key falls after the last bucket
      {
          bucket = spline_points_.getBucket(*(--it));
          rs::BucketEntry<KeyType>* left = bucket->getLastElement(); //first entry of the right bucket
          if (left->getRightBufferSize() == overflow_buffer_size)
            return retrain();
          else
            left->insertRight(newKey, pos);
      }
    }

    else
    {
      size_t entry_pos = bucket->get_entry_pos_equal_or_lowest_greater_than_key(newKey);
      size_t bucket_size = bucket->getNumElements();
      rs::BucketEntry<KeyType>* entry = bucket->getElement(entry_pos);
      if(entry_pos > 0 && entry_pos < bucket_size) //if newKey belongs to this bucket
      {
        if (entry->getLeftBufferSize() == overflow_buffer_size)
          return retrain();
        else  
          entry->insertLeft(newKey, pos);
      }
      else if (entry_pos == 0)  //if newKey is less than the first entry of this bucket
      {
        if (entry->getLeftBufferSize() < overflow_buffer_size)
          entry->insertLeft(newKey, pos);
        else
        {
          auto it = std::lower_bound(radix_table_.begin(), radix_table_.end(), prefix);
          if (it == radix_table_.begin())
          {
            return retrain();
          }
          else
          {
            bucket = spline_points_.getBucket(*(--it));
            rs::BucketEntry<KeyType>* entry = bucket->getLastElement();
            if(entry->getRightBufferSize() < overflow_buffer_size)
              entry->insertRight(newKey, pos);
            else
              return retrain();  
          }
        }
      }
      else   //if newKey is greater than the last entry of this bucket
      {
        entry = bucket->getLastElement();
        if (entry->getRightBufferSize() < overflow_buffer_size)
          entry->insertRight(newKey, pos);
        else
        {
          auto it = std::lower_bound(radix_table_.begin(), radix_table_.end(), prefix);
          if (it == radix_table_.end())
          {
            return retrain();
          }
          else
          {
            bucket = spline_points_.getBucket(*(++it));
            rs::BucketEntry<KeyType>* entry = bucket->getElement(0);
            if(entry->getLeftBufferSize() < overflow_buffer_size)
              entry->insertLeft(newKey, pos);
            else
              return retrain();  
          }
        }
      }
    }
    return 1;

  }


  std::pair<Coord<KeyType>*, size_t> getAllInsertedData()
  {
    rs::Coord<KeyType>* resInBuffer = nullptr;
    size_t resSize = 0;
    rs::Bucket<KeyType> *curr_bucket = nullptr, *prev_bucket = nullptr;

    for(int i = 0 ; i < radix_table_.size() ; i++)
    {
      KeyType prefix = radix_table_[i];
      curr_bucket = spline_points_.getBucket(prefix);
      if (curr_bucket == prev_bucket)
        continue;
      
      for (int i = 0 ; i < curr_bucket->getNumElements() ; i++)
      {
        rs::BucketEntry<KeyType>* entry = curr_bucket->getElement(i);
        size_t tempSize = entry->getLeftBufferSize();
        if (tempSize > 0)
        {
          
          rs::Coord<KeyType>* temp = new rs::Coord<KeyType> [resSize + tempSize];
          std::copy(resInBuffer, resInBuffer + resSize, temp);
          std::copy(entry->left_overflow_buffer, entry->left_overflow_buffer+tempSize, temp+resSize);
          resSize += tempSize;
          resInBuffer = temp;
        }

        tempSize = entry->getRightBufferSize();
        if (tempSize > 0)
        {
          
          rs::Coord<KeyType>* temp = new rs::Coord<KeyType> [resSize + tempSize];
          std::copy(resInBuffer, resInBuffer + resSize, temp);
          std::copy(entry->right_overflow_buffer, entry->right_overflow_buffer+tempSize, temp+resSize);
          resSize += tempSize;
          resInBuffer = temp;
        }

        prev_bucket = curr_bucket;
      }

    }

    return {resInBuffer, resSize};
  }


 private:
  // Returns the index of the spline point that marks the end of the spline
  // segment that contains the `key`: `key` ∈ (spline[index - 1], spline[index]]

  /*
  size_t GetSplineSegment(const KeyType key) const {
    // Narrow search range using radix table.
    const KeyType prefix = (key - min_key_) >> num_shift_bits_;
    assert(prefix + 1 < radix_table_.size());
    const uint32_t begin = radix_table_[prefix];
    const uint32_t end = radix_table_[prefix + 1];

    if (end - begin < 32) {
      // Do linear search over narrowed range.
      uint32_t current = begin;
      while (spline_points_[current].x < key) ++current;
      return current;
    }

    // Do binary search over narrowed range.
    const auto lb = std::lower_bound(
        spline_points_.begin() + begin, spline_points_.begin() + end, key,
        [](const Coord<KeyType>& coord, const KeyType key) {
          return coord.x < key;
        });
    return std::distance(spline_points_.begin(), lb);
  }
  */

  static bool compare_by_key (const Coord<KeyType>& p1, const Coord<KeyType>& p2) {
    return p1.x < p2.x;
  }
  

  std::pair<rs::BucketEntry<KeyType>*, rs::BucketEntry<KeyType>*> GetSplineSegment(const KeyType key) const {
    // Narrow search range using radix table.
    //const KeyType prefix = (key - min_key_) >> num_shift_bits_;
    const KeyType prefix = spline_points_.getPrefix(key);
    auto it = std::lower_bound(radix_table_.begin(), radix_table_.end(), prefix);
    KeyType bucket_id1, bucket_id2;
    rs::BucketEntry<KeyType> *firstEntry, *secondEntry;

    bool flag = false;

    if (*it == prefix) //if there is a bucket for the prefix
    {
      rs::Bucket<KeyType> *temp_bucket = spline_points_.getBucket(*it);
      rs::BucketEntry<KeyType> *minEntry = temp_bucket->getElement(0);
      rs::BucketEntry<KeyType> *maxEntry = temp_bucket->getLastElement();
      if (key < minEntry->val) //if key lies in between this bucket and previous bucket
      {
        secondEntry = minEntry;
        temp_bucket = spline_points_.getBucket(*(--it));
        rs::BucketEntry<KeyType>* p = temp_bucket->getLastElement();
        firstEntry = p;
        return std::pair<rs::BucketEntry<KeyType>*, rs::BucketEntry<KeyType>*>(firstEntry, secondEntry);    
      }
      else if (key > maxEntry->val) //if key lies in between this bucket and next bucket
      {
        firstEntry = maxEntry;
        temp_bucket = spline_points_.getBucket(*(++it));
        rs::BucketEntry<KeyType>* p = temp_bucket->getElement(0);
        secondEntry = p;
        return std::pair<rs::BucketEntry<KeyType>*, rs::BucketEntry<KeyType>*>(firstEntry, secondEntry);    
      }
      else //key falls within this bucket
      {
        bucket_id1 = bucket_id2 = prefix;
        flag == true;
      }
    }
    else
    {
      it--;
      bucket_id1 = *it;
      bucket_id2 = *(it+1);
    } 
    
    //if key falls in the same bucket
    if (bucket_id1 == bucket_id2)
    {
      rs::Bucket<KeyType> *temp_bucket = spline_points_.getBucket(bucket_id1);
      rs::BucketEntry<KeyType> target(key, -1);
      auto p = temp_bucket->get_equal_or_highest_less_than_key(key);
      firstEntry = &(*p);

      if(firstEntry->val == key)
      {
        secondEntry = &(*p);
      }
      else
      {
        p++;
        secondEntry = &(*p);
      }     
    }
    else
    {
      rs::Bucket<KeyType> *temp_bucket = spline_points_.getBucket(bucket_id1);
      rs::BucketEntry<KeyType>* p = temp_bucket->getLastElement();
      firstEntry = p;

      temp_bucket = spline_points_.getBucket(bucket_id2);
      p = temp_bucket->getElement(0);
      secondEntry = p;
    }

    return std::pair<rs::BucketEntry<KeyType>*, rs::BucketEntry<KeyType>*>(firstEntry, secondEntry);    

  }

  int retrain()
  {
    //std::cout << "Model Retraining required" << std::endl;
    return -1;
  }


  KeyType min_key_;
  KeyType max_key_;
  size_t num_keys_;
  size_t num_radix_bits_;
  size_t num_shift_bits_;
  size_t max_error_;
  size_t overflow_buffer_size;

  std::vector<uint32_t> radix_table_;
  //std::vector<rs::Coord<KeyType>> spline_points_;   //skd
  rs::Directory<KeyType>spline_points_;

  template <typename>
  friend class Serializer;
};

}  // namespace rs