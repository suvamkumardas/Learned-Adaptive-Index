#pragma once

#include <bits/stdc++.h>
#include <iterator>
namespace rs {

template <typename KeyType>
class BucketEntry {
  public:
    KeyType val;
    double pos;
    rs::Coord<KeyType> *left_overflow_buffer, *right_overflow_buffer;
    size_t left_size, right_size, max_buffer_size;

     BucketEntry() {}

     BucketEntry(KeyType val, double pos, Coord<KeyType> *left, Coord<KeyType> *right, 
        size_t max_buffer_size, size_t left_buffer_size, size_t right_buffer_size)
     {
        this->val = val;
        this->pos = pos;
        this->left_overflow_buffer = left;
        this->right_overflow_buffer = right;
        this->max_buffer_size = max_buffer_size;
        this->left_size = left_buffer_size;
        this->right_size = right_buffer_size;
     }
     BucketEntry(KeyType val, double pos, int max_buffer_size) : BucketEntry(val, pos, nullptr, nullptr, max_buffer_size, 0, 0)
     {
     }

     BucketEntry(KeyType val, double pos) : BucketEntry(val, pos, 0)
     {
     }

     void insertLeft(KeyType val, double pos)
     {
        assert(left_size < max_buffer_size);
        if (left_size == 0)
            left_overflow_buffer = new rs::Coord<KeyType> [max_buffer_size];
        left_overflow_buffer[left_size] = {val, pos};
        left_size++;   
     }

     void insertRight(KeyType val, double pos)
     {
        assert(right_size < max_buffer_size);
        if (right_size == 0)
            right_overflow_buffer = new rs::Coord<KeyType> [max_buffer_size];
        right_overflow_buffer[right_size] = {val, pos};
        right_size++;           
     }

     size_t getLeftBufferSize()
     {
        return left_size;
     }

     size_t getRightBufferSize()
     {
        return right_size;
     }
};

template <typename KeyType>
class Bucket {
        int depth,size, max_size, overflow_buffer_size;
        BucketEntry<KeyType>* values;

#if 0
    //Custom comparator to compare by key
    bool compare_by_key = [](const std::pair<KeyType, double>& p1, const std::pair<KeyType, double>& p2) {
        return p1.first < p2.first;
    };
#endif    

    //Custom comparator to compare by key
    static bool compare_by_key (const BucketEntry<KeyType>& p1, const BucketEntry<KeyType>& p2) {
        return p1.val < p2.val;
    };

    size_t find (KeyType key)
    {
        BucketEntry<KeyType> target(key, -1);
        auto it = std::lower_bound(values, values+size, target, compare_by_key);
        
        size_t pos = it - values;
        if (pos >= size)
            return -1;
        else if (values[pos].val == key)
            return pos;
        else                    //if (values[pos].first != key)
            return -1;     
    }


    public:

        Bucket(int depth, int size, int overflow_buffer_size)
        {
            this->depth = depth;
            this->max_size = size;
            this->size = 0;
            values = new BucketEntry<KeyType> [max_size];
            this->overflow_buffer_size = overflow_buffer_size;
        }

        int insert(KeyType key, double value)
        {
            //std::map<KeyType,double>::iterator it;
            size_t pos = find(key);
            if(pos != -1)
                return -1;
            if(isFull())
                return 0;
            values[size] = BucketEntry(key, value, overflow_buffer_size);
            size++;
            return 1;
        }

#if 0
        int remove(KeyType key)
        {
            //std::map<KeyType,double>::iterator it;
            auto it = values.find(key);
            if(it!=values.end())
            {
                values.erase(it);
                return 1;
            }
            else
            {
                std::cout<<"Cannot remove : This key does not exists"<<std::endl;
                return 0;
            }
        }
        

        int update(KeyType key, double value)
        {
            //std::map<KeyType,double>::iterator it;
            auto it = values.find(key);
            if(it!=values.end())
            {
                values[key] = value;
                std::cout<<"Value updated"<<std::endl;
                return 1;
            }
            else
            {
                std::cout<<"Cannot update : This key does not exists"<<std::endl;
                return 0;
            }
        }
#endif        

        size_t search(KeyType key)
        {
            //std::map<KeyType,KeyType>::iterator it;
            size_t pos = find(key);
            return pos;
        }

        int isFull(void)
        {
            if(size == max_size)
                return 1;
            else
                return 0;
        }

        int isEmpty(void)
        {
            if(size==0)
                return 1;
            else
                return 0;
        }

        int getDepth(void)
        {
            return depth;
        }

        int increaseDepth(void)
        {
            depth++;
            return depth;
        }

        int decreaseDepth(void)
        {
            depth--;
            return depth;
        }

        BucketEntry<KeyType>* copy(void)
        {
            BucketEntry<KeyType>* temp = new BucketEntry<KeyType> [max_size];
            std::copy(values, values+max_size, temp);
            return temp;
        }

        void clear(void)
        {
            //std::fill(values, values + max_size, std::make_pair(0, -1.0));
            //size = 0;
            if (size > 0)
            {
                delete [] values;
                size = 0;
            } 
        }

        void display()
        {
            //std::map<KeyType,KeyType>::iterator it;
            for(auto it=std::begin(values);it!=std::end(values);it++)
                std::cout<<it->val<<" ";
            std::cout<<std::endl;
        }

        int getNumElements()
        {
            return size;
        }

#if 0
        std::pair<KeyType, double> getLastKeyValuePair()
        {
            auto lastPair = values.end();
            //lastPair--;
            return {lastPair->first, lastPair->second};
        }

        std::pair<KeyType, double>* getMap(void)
        {
            return values;
        }
#endif        

        auto get_equal_or_highest_less_than_key(KeyType key)
        {
            BucketEntry<KeyType> target(key, -1);
            auto it = std::lower_bound(values, values+size, target, compare_by_key);
            return (*it).val == key ? it : it-1;
        }

        size_t get_entry_pos_equal_or_lowest_greater_than_key(KeyType key)
        {
            BucketEntry<KeyType> target(key, -1);
            auto it = std::lower_bound(values, values+size, target, compare_by_key);
            size_t pos = it - values;
            return (pos < size) ? pos : -1;
        }

        BucketEntry<KeyType>* end()
        {
            return std::end(values);
        }

        BucketEntry<KeyType>* getElement(int pos)
        {
            return &values[pos];
        }

        BucketEntry<KeyType>* getLastElement()
        {
            return &values[size-1];
        }

};

template <typename KeyType>
class Directory {
        int initial_depth, global_depth, bucket_size, dir_size, num_bits_to_shift, overflow_buffer_size;
        //std::vector<Bucket<KeyType>*> buckets;
        Bucket<KeyType>** buckets;
        KeyType min_key;
        
        KeyType hash(KeyType n) const
        {
            //return n&((1<<global_depth)-1);
            //const uint32_t clzl = __builtin_clzl(n - min_key);
            //size_t numBitsToShift =  sizeof(n) * 8 - clzl - global_depth;
            //return (n - min_key) >> (sizeof(n) * 8 - global_depth);
            
            KeyType temp = (n - min_key) >> num_bits_to_shift;
            //return temp&((1<<global_depth)-1);
            return temp;

            #if 0
            KeyType temp = (n - min_key) >> num_bits_to_shift;
            const uint32_t to_shift_mask = sizeof(n) * 8 - global_depth - num_bits_to_shift;
            KeyType mask = (1<<global_depth)-1;
            mask = mask << to_shift_mask;
            return temp & mask;
            #endif
        }

        int pairIndex(int bucket_no, int depth)
        {
            return bucket_no^(1<<(depth-1));
        }

#if 0
        void grow(void)
        {
            for(int i = 0 ; i < 1<<global_depth ; i++ )
                buckets.push_back(buckets[i]);
            global_depth++;
        }
#endif        
        void grow(void)
        {
            int new_dir_size = dir_size << 1;
            Bucket<KeyType>** new_buckets = new Bucket<KeyType>* [new_dir_size];
            for (int i = 0 ; i < dir_size ; i++)
            {
                int idx = i << 1;
                new_buckets[idx] = buckets[i]; 
                new_buckets[idx+1] = buckets[i];           
            }
            delete[] buckets;
            buckets = new_buckets;
            global_depth++;
            num_bits_to_shift--;
            dir_size = new_dir_size;
            std::cout << "Global Depth Increased = " << global_depth << std::endl;
        }


#if 0
        void shrink(void)
        {
            int i,flag=1;
            for( i=0 ; i<buckets.size() ; i++ )
            {
                if(buckets[i]->getDepth()==global_depth)
                {
                    flag=0;
                    return;
                }
            }
            global_depth--;
            for(i = 0 ; i < 1<<global_depth ; i++ )
                buckets.pop_back();
        }
#endif        

#if 0
        void split(int bucket_no)
        {
            int local_depth,pair_index,index_diff,dir_size,i;
            std::pair<KeyType, double>* temp;
            //map<KeyType, double>::iterator it;

            local_depth = buckets[bucket_no]->increaseDepth();
            if(local_depth>global_depth)
                grow();
            pair_index = pairIndex(bucket_no,local_depth);
            buckets[pair_index] = new Bucket<KeyType>(local_depth,bucket_size);
            temp = buckets[bucket_no]->copy();
            int temp_bucket_size = buckets[bucket_no]->getNumElements();
            buckets[bucket_no]->clear();
            index_diff = 1<<local_depth;
            dir_size = 1<<global_depth;
            for( i=pair_index-index_diff ; i>=0 ; i-=index_diff )
                buckets[i] = buckets[pair_index];
            for( i=pair_index+index_diff ; i<dir_size ; i+=index_diff )
                buckets[i] = buckets[pair_index];
            for(int it=0 ; it <= temp_bucket_size ; it++)
                insert((temp[it]).first,(temp[it]).second,1);
        }
#endif        

        void reinitialize_bucket(int bucket_num, int local_depth)
        {
            //buckets[bucket_num]->clear();
            //delete buckets[bucket_num];
            //buckets[bucket_num] = nullptr;
            buckets[bucket_num] = new Bucket<KeyType>(local_depth,bucket_size,overflow_buffer_size);
        }

        void split(int bucket_no)
        {
            //map<KeyType, double>::iterator it;

            int local_depth = buckets[bucket_no]->increaseDepth();
            BucketEntry<KeyType>* temp = buckets[bucket_no]->copy();
            int temp_bucket_size = buckets[bucket_no]->getNumElements();
            int index1, index2;
            if(local_depth>global_depth)
            {
                grow();
                index1 = bucket_no << 1;
                index2 = index1 + 1;
                reinitialize_bucket(index1, local_depth);
                reinitialize_bucket(index2, local_depth);
            }
            else
            {
                if (!(bucket_no & 1))
                {
                    index1 = bucket_no;
                    index2 = bucket_no + 1;
                }
                else
                {
                    index1 = bucket_no - 1;
                    index2 = bucket_no;
                }

                if (buckets[index1] != buckets[index2])
                    reinitialize_bucket(bucket_no, local_depth);
                else
                {
                    reinitialize_bucket(index1, local_depth);
                    reinitialize_bucket(index2, local_depth);
                }          
            }
            //reinitialize_bucket(index1, local_depth);
            //reinitialize_bucket(index2, local_depth);
            
            for(int it=0 ; it < temp_bucket_size ; it++)
                insert((temp[it]).val,(temp[it]).pos,1);
        }


        void merge(int bucket_no)
        {
            int local_depth,pair_index,index_diff,dir_size,i;

            local_depth = buckets[bucket_no]->getDepth();
            pair_index = pairIndex(bucket_no,local_depth);
            index_diff = 1<<local_depth;
            dir_size = 1<<global_depth;

            if( buckets[pair_index]->getDepth() == local_depth )
            {
                buckets[pair_index]->decreaseDepth();
                delete(buckets[bucket_no]);
                buckets[bucket_no] = buckets[pair_index];
                for( i=bucket_no-index_diff ; i>=0 ; i-=index_diff )
                    buckets[i] = buckets[pair_index];
                for( i=bucket_no+index_diff ; i<dir_size ; i+=index_diff )
                    buckets[i] = buckets[pair_index];
            }
        }

        std::string bucket_id(int n, int global_depth)
        {
            int d;
            std::string s;
            //d = buckets[n]->getDepth();
            d = global_depth;
            s = "";
            while(n>0 && d>0)
            {
                s = (n%2==0?"0":"1")+s;
                n/=2;
                d--;
            }
            while(d>0)
            {
                s = "0"+s;
                d--;
            }
            return s;
        }

    public:

        Directory() {}
        
        Directory(KeyType min_key, int depth, int num_bits_to_shift, int bucket_size, int overflow_buffer_size)
        {
            this->min_key = min_key;
            this->global_depth = this->initial_depth = depth;
            this->num_bits_to_shift = num_bits_to_shift;
            this->bucket_size = bucket_size;
            this->dir_size = 1<<depth;
            this->overflow_buffer_size = overflow_buffer_size;
            buckets = new Bucket<KeyType>* [this->dir_size];
            std::fill(buckets, buckets+this->dir_size, nullptr);
        #if 0    
            for(int i = 0 ; i < this->dir_size ; i++ )
            {
                buckets[i] = new Bucket<KeyType>(depth,bucket_size);
            }
        #endif    
        }

        KeyType insert(KeyType key,double value,bool reinserted)
        {            
            KeyType bucket_no = hash(key);
            
            if (buckets[bucket_no] == nullptr)
                buckets[bucket_no] = new Bucket<KeyType>(initial_depth,bucket_size, overflow_buffer_size);
            
            //Bucket<KeyType>* ptr_array[256];
            //for int(i = 0 ; i < buckets[bucket_no];
            int status = buckets[bucket_no]->insert(key,value);
#if 0
            if(status==1)
            {
                if(!reinserted)
                    std::cout<<"Inserted key "<<key<<" in bucket "<<bucket_id(bucket_no, global_depth) << " - " << bucket_no <<std::endl;
                else
                    std::cout<<"Moved key "<<key<<" to bucket "<<bucket_id(bucket_no, global_depth) << " - " << bucket_no <<std::endl;
            }
            else if(status==0)
            {
                split(bucket_no);
                bucket_no = insert(key,value,reinserted);
            }
            else
            {
                std::cout<<"Key "<<key<<" already exists in bucket "<<bucket_id(bucket_no, global_depth)<<std::endl;
            }
#else            

            if(status==0)
            {
                split(bucket_no);
                bucket_no = insert(key,value,reinserted);
            }
#endif            
            return bucket_no;
        }

#if 0
        void remove(KeyType key,int mode)
        {
            int bucket_no = hash(key);
            if(buckets[bucket_no]->remove(key))
                std::cout<<"Deleted key "<<key<<" from bucket "<<bucket_id(bucket_no)<<std::endl;
            if(mode>0)
            {
                if(buckets[bucket_no]->isEmpty() && buckets[bucket_no]->getDepth()>1)
                    merge(bucket_no);
            }
            if(mode>1)
            {
                shrink();
            }
        }
        

        void update(KeyType key, double value)
        {
            int bucket_no = hash(key);
            buckets[bucket_no]->update(key,value);
        }
#endif

        size_t search(KeyType key)
        {
            int bucket_no = hash(key);
            //std::cout<<"Searching key "<<key<<" in bucket "<<bucket_id(bucket_no)<<std::endl;
            size_t pos = buckets[bucket_no]->search(key);
            return pos;
        }

        void display(bool duplicates)
        {
            int i,j,d;
            std::string s;
            std::set<std::string> shown;
            std::cout<<"Global depth : "<<global_depth<<std::endl;
            for(i=0;i<dir_size;i++)
            {
                d = buckets[i]->getDepth();
                s = bucket_id(i);
                if(duplicates || shown.find(s)==shown.end())
                {
                    shown.insert(s);
                    for(j=d;j<=global_depth;j++)
                        std::cout<<" ";
                    std::cout<<s<<" => ";
                    buckets[i]->display();
                }
            }
        }

        int globalElementCount()
        {
            int size = 0;
            for (int i = 0 ; i < dir_size ; i++)
            {
                size += buckets[i]->getNumElements();
            }
            return size;
        }

#if 0
        std::pair<KeyType, double> getLastKeyValuePair(void)
        {
            Bucket<KeyType> *last = buckets[dir_size - 1];
            return last->getLastKeyValuePair();
        }
#endif        

        int getDirectoryLength() const
        {
            return this->dir_size;
        }

        Bucket<KeyType>* getBucket(KeyType pos) const
        {
            return buckets[pos];
        }

#if 0
        int getBucketNumber(KeyType key) const
        {
            int bucket_no = key &((1<<global_depth)-1);
            return bucket_no;
        }
#endif        

        KeyType getPrefix(KeyType n) const
        {
            int prefix = hash(n);
            return prefix;
        }

        std::vector<unsigned int> getAllPrefixes()
        {
            std::vector<unsigned int> v;
            v.reserve(dir_size);
            for (unsigned int i = 0 ; i < dir_size ; i++)
                if (buckets[i] != nullptr && buckets[i]->getNumElements() > 0)
                    v.push_back(i);

            return v;        
        }

        size_t getOverflowBufferSize()
        {
            return overflow_buffer_size;
        }

};


//void menu();

/* Main function */
/*

int main()
{
    bool show_messages, show_duplicate_buckets;
    int bucket_size, initial_global_depth, value;
    int key, mode;
    string choice;

    // Set show_messages to 0 when taking input using file redirection
    show_messages = 0;

    // Set show_duplicate_buckets to 1 to see all pointers instead of unique ones
    show_duplicate_buckets = 0;

    if(show_messages) { cout<<"Bucket size : "; }
    cin>>bucket_size;
    if(show_messages) { cout<<"Initial global depth : "; }
    cin>>initial_global_depth;

    Directory<int> d(initial_global_depth,bucket_size);
    cout<<endl<<"Initialized directory structure"<<endl;

    if(show_messages)
        menu();

    do
    {
        cout<<endl;
        if(show_messages) { cout<<">>> "; }
        cin>>choice;
        if(choice=="insert")
        {
            cin>>key>>value;
            if(show_messages) { cout<<endl; }
            d.insert(key,value,0);
        }
        else if(choice=="delete")
        {
            cin>>key>>mode;
            if(show_messages) { cout<<endl; }
            d.remove(key,mode);
        }
        else if(choice=="update")
        {
            cin>>key>>value;
            if(show_messages) { cout<<endl; }
            d.update(key,value);
        }
        else if(choice=="search")
        {
            cin>>key;
            if(show_messages) { cout<<endl; }
            d.search(key);
        }
        else if(choice=="display")
        {
            if(show_messages) { cout<<endl; }
            d.display(show_duplicate_buckets);
        }
    } while(choice!="exit");

    return 0;
}
*/

/* Print usage menu */
/*
void menu()
{
    cout<<"--------------------"<<endl;
    cout<<"Enter queries in the following format :"<<endl;
    cout<<"insert <key> <value>     (key: integer, value: string)"<<endl;
    cout<<"delete <key> <mode>      (mode: 0-Lazy / 1-Merge empty buckets / 2-Merge buckets and shrink )"<<endl;
    cout<<"update <key> <new value>"<<endl;
    cout<<"search <key>"<<endl;
    cout<<"display"<<endl;
    cout<<"exit"<<endl;
    cout<<"--------------------"<<endl;
}
*/

}
