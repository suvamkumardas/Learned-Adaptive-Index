#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <bits/stdc++.h>
#include <chrono>
#include <unistd.h>
#include <filesystem>

#include "LearnedAdaptiveIndex/LAI.hpp"


int main(int argc, char* argv[])
{
    if (argc < 4)
    {
        std::cerr << "Usage : LAI <data_file> <dataset_size> <query_workload> <sub-runs> [max-error-bound]" << std::endl;
        exit(0);
    }
    
    using datatype_t = unsigned long;
    std::string base_path = std::string(std::filesystem::current_path()) + "/data_and_queries/100_million/";
    std::cout << "Base Path = " << base_path << std::endl;
    std::ifstream data_file(base_path + argv[1]);
    int inputDataSize = atoi(argv[2]), subRuns = atoi(argv[4]);
    size_t error_bound = argc == 6 ? atoi(argv[5]) : 32;
    datatype_t* dataset =  new datatype_t [inputDataSize + 2];
    datatype_t key;

    int pos = 1;
    while(data_file >> key)
    {
        dataset[pos] = key;
        pos++;
    }
    std::cout << "Dataset of size: " << inputDataSize << std::endl;
    std::cout << "Pos: " << pos << std::endl;
    assert(pos - 1 == inputDataSize);
    std::cout << "Running LAI without predictions" << std::endl;
    std::cout << "Learned Index Error Bound: " << error_bound << std::endl;
    std::cout.flush();

    LAI <datatype_t> lai(dataset, inputDataSize, error_bound, false);
    double totalTimeElapsed = 0, queryTimeElapsed = 0, predictionTimeElapsed = 0, cumulativeQueryTimeOnly = 0;

    std::string query_folder = base_path + "queries_sel_random/";
    std::string real_query_file_path = query_folder + argv[3] + ".txt";
    std::ifstream query_file(real_query_file_path);

    if (!query_file.is_open())
    {
        std::cout << "Query file is missing.. Abort.." <<endl;
        exit(0);
    }

    std::cout << "Query File Path = " << real_query_file_path << endl;
    std::string line;
    int query_num = 0;
    std::cout << "\nquery_num,Cumulative Query Time Only,Total Time Elapsed" << endl; 
    while (query_file >> line)
    {
        char *query_bounds = line.data();
        char *token = strtok(query_bounds, ",");
        datatype_t lowVal = atoi(token);
        token = strtok(NULL, ",");
        datatype_t highVal = atoi(token);
        
        auto started = std::chrono::high_resolution_clock::now();    
        std::pair<long,long> pos = lai.query(lowVal,highVal);
        auto done = std::chrono::high_resolution_clock::now();
        auto lowValPos = pos.first, highValPos = pos.second;
        double time = (double)std::chrono::duration_cast<std::chrono::nanoseconds>(done-started).count();
        
        queryTimeElapsed += time;
        query_num++;
        if(query_num == 1 || (query_num % subRuns == 0))
        {
            std::cout << queryTimeElapsed/pow(10,9) << "," << queryTimeElapsed/pow(10,9) << endl;
            std::cout.flush();
        }
}
    lai.getQueryStats();
}

