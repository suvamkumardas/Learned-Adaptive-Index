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
        spline_points_(std::move(spline_points))
        {
        }

  // Returns the estimated position of `key`.
  double GetEstimatedPosition(const KeyType key) const {
    // Truncate to data boundaries.
    if (key <= min_key_) return 0;
    if (key >= max_key_) return num_keys_ - 1;

    
    // Find spline segment with `key` ∈ (spline[index - 1], spline[index]].
    /*  //skd
    const size_t index = GetSplineSegment(key);
    const Coord<KeyType> down = spline_points_[index - 1];
    const Coord<KeyType> up = spline_points_[index];
    */
    const std::pair<Coord<KeyType>, Coord<KeyType>> p = GetSplineSegment(key);
    const Coord<KeyType> down = p.first;
    const Coord<KeyType> up = p.second;

    // Compute slope.
    const double x_diff = up.x - down.x;
    const double y_diff = up.y - down.y;
    const double slope = y_diff / x_diff;

    // Interpolate.
    const double key_diff = key - down.x;
    return std::fma(key_diff, slope, down.y);
  }

  // Returns a search bound [begin, end) around the estimated position.
  SearchBound GetSearchBound(const KeyType key) const {
    const size_t estimate = GetEstimatedPosition(key);
    const size_t begin = (estimate < max_error_) ? 0 : (estimate - max_error_);
    // `end` is exclusive.
    const size_t end = (estimate + max_error_ + 2 > num_keys_)
                           ? num_keys_
                           : (estimate + max_error_ + 2);
    return SearchBound{begin, end};
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

  std::pair<Coord<KeyType>,Coord<KeyType>> GetSplineSegment(const KeyType key) const {
    // Narrow search range using radix table.
    //const KeyType prefix = (key - min_key_) >> num_shift_bits_;
    const KeyType prefix = spline_points_.getPrefix(key);
    auto it = std::lower_bound(radix_table_.begin(), radix_table_.end(), prefix);
    KeyType bucket_id1, bucket_id2;
    KeyType firstKey, secondKey;
    double firstVal, secondVal;
    bool flag = false;

    if (*it == prefix) //if there is a bucket for the prefix
    {
      rs::Bucket<KeyType> *temp_bucket = spline_points_.getBucket(*it);
      rs::BucketEntry<KeyType> *minEntry = temp_bucket->getElement(0);
      rs::BucketEntry<KeyType> *maxEntry = temp_bucket->getLastElement();
      if (key < minEntry->val) //if key lies in between this bucket and previous bucket
      {
        secondKey = minEntry->val;
        secondVal = minEntry->pos;
        temp_bucket = spline_points_.getBucket(*(--it));
        rs::BucketEntry<KeyType>* p = temp_bucket->getLastElement();
        firstKey = p->val;
        firstVal = p->pos;
        Coord<KeyType> c1 = {firstKey, firstVal}, c2 = {secondKey, secondVal};
        return std::pair<Coord<KeyType>, Coord<KeyType>>(c1, c2);    
      }
      else if (key > maxEntry->val) //if key lies in between this bucket and next bucket
      {
        firstKey = maxEntry->val;
        firstVal = maxEntry->pos;
        temp_bucket = spline_points_.getBucket(*(++it));
        rs::BucketEntry<KeyType>* p = temp_bucket->getElement(0);
        secondKey = p->val;
        secondVal = p->pos;
        Coord<KeyType> c1 = {firstKey, firstVal}, c2 = {secondKey, secondVal};
        return std::pair<Coord<KeyType>, Coord<KeyType>>(c1, c2);    
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
      firstKey = (*p).val;
      firstVal = (*p).pos;

      if(firstKey == key)
      {
        secondKey = p->val;
        secondVal = p->pos;
      }
      else
      {
        p++;
        secondKey = (*p).val;
        secondVal = (*p).pos;
      }     
    }
    else
    {
      rs::Bucket<KeyType> *temp_bucket = spline_points_.getBucket(bucket_id1);
      rs::BucketEntry<KeyType>* p = temp_bucket->getLastElement();
      firstKey = p->val;
      firstVal = p->pos;

      temp_bucket = spline_points_.getBucket(bucket_id2);
      p = temp_bucket->getElement(0);
      secondKey = p->val;
      secondVal = p->pos;
    }

    Coord<KeyType> c1 = {firstKey, firstVal}, c2 = {secondKey, secondVal};
    return std::pair<Coord<KeyType>, Coord<KeyType>>(c1, c2);    

#if 0    
    rs::Bucket<KeyType> *temp_bucket = spline_points_.getBucket(bucket_id1);    

    std::pair<KeyType, double>* last_element = temp_bucket->getLastElement();
    std::pair<KeyType, double>* lb = flag ? temp_bucket->get_highest_less_than_key(key)
                                          : last_element;
    //auto temp = lb;
    //lb--;
    if (lb->first == key)
    {
        firstKey = secondKey = key;
        firstVal = secondVal = lb->second;        
        Coord<KeyType> c1 = {firstKey, firstVal}, c2 = {secondKey, secondVal};
        return std::pair<Coord<KeyType>, Coord<KeyType>>(c1, c2);    
    }
    
    firstKey = lb->first;  
    firstVal = lb->second;

    if (!flag && lb != last_element)
    {
        secondKey = (lb+1)->first;  
        secondVal = (lb+1)->second;
    }
    else
    {
        rs::Bucket<KeyType> *temp_bucket2 = spline_points_.getBucket(bucket_id2);
        std::pair<KeyType, double> p = temp_bucket2->getElement(0);
        secondKey = p.first;  
        secondVal = p.second;
    }
        
    Coord<KeyType> c1 = {firstKey, firstVal}, c2 = {secondKey, secondVal};
    return std::pair<Coord<KeyType>, Coord<KeyType>>(c1, c2);
#endif

  }

  int retrain()
  {
    std::cout << "Model Retraining required" << std::endl;
    return -1;
  }

  int insertNewKey(KeyType newKey, double pos)
  {
    const KeyType prefix = spline_points_.getPrefix(newKey);
    rs::Bucket<KeyType>* bucket = spline_points_.getBucket(prefix);
    
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
      if(entry_pos > 0 && entry_pos < bucket_size-1) //if newKey belongs to this bucket
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
        assert(entry_pos == bucket_size-1);
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
            rs::BucketEntry<KeyType>* entry = bucket->Element(0);
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