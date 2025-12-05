#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <bits/stdc++.h>
#include <random>
#include "LearnedAdaptiveIndex/LAI.hpp" 

int main (int argc, char *argv[])
{
    using datatype_t = unsigned int;
    std::ifstream data_file("data_and_queries/YCSB2/input_dataset_100M_unique.txt");
    std::ifstream query_file(argv[1]);
    int subRuns = 2000;
    //std::ofstream output_file("data_and_queries/YCSB/input_dataset_100M_unique.txt");
    std::vector<datatype_t> dataset;
    datatype_t key, low, high, insertkey;

    while(data_file >> key)
        dataset.push_back(key);
    
    int dataset_size =  dataset.size();
    std::cout << "Initial Dataset Size = " << dataset_size << "\n";
    
    std::vector<rs::Coord<datatype_t>> coord_vec;
    for(unsigned int i = 0 ; i < dataset.size() ; i++)
    {
        rs::Coord<datatype_t> c = {dataset[i], i};
        coord_vec.push_back(c);
    }
   
    LAI <datatype_t> lai(coord_vec, 18, 32, 128, true);
    
    std::cout << query_file.is_open() << " Path = " << argv[1] << std::endl;
    std::string line;
    int query_num = 0;
    double totalTimeElapsed = 0, time;

    while (query_file >> line)
    {
        int size = 0;
        if (line.find(",") != std::string::npos)
        {
            char *query_bounds = line.data();
            char *token = strtok(query_bounds, ",");
            datatype_t lowVal = atoi(token);
            token = strtok(NULL, ",");
            datatype_t highVal = atoi(token);

            /*
            for(int i = 0 ; i < dataset.size() ; i++)
                {
                    if (lowVal <= dataset[i] && dataset[i] <= highVal)
                        size++;
                }
            */            

            auto started = std::chrono::high_resolution_clock::now();    
            std::tuple<long,long, rs::Coord<unsigned int>*, size_t> pos = lai.query(lowVal,highVal);
            auto done = std::chrono::high_resolution_clock::now();
            time = (double)std::chrono::duration_cast<std::chrono::nanoseconds>(done-started).count();           
        }
        else
        {
            datatype_t insertKey = std::stoul(line);
            rs::Coord<datatype_t> c = {insertKey, 0};
            std::vector<rs::Coord<datatype_t>> v = {c};
            auto started = std::chrono::high_resolution_clock::now();    
            lai.updatePendingInserts(v);
            auto done = std::chrono::high_resolution_clock::now();
            time = (double)std::chrono::duration_cast<std::chrono::nanoseconds>(done-started).count();           
        }
        totalTimeElapsed += time;
        query_num++;
    
        if(query_num == 1 || query_num % subRuns == 0)
        {
            std::cout << "Cumulative Runtime (" << query_num << " queries)," << totalTimeElapsed/pow(10,9) << endl; //in seconds
            std::cout.flush();
            //std::cout << "Res Size = " << size;
            //break;
        }
    }
    lai.getQueryStats();
}
