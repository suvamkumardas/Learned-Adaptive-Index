#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <bits/stdc++.h>
#include <chrono>
#include <unistd.h>
#include <filesystem>

#include "LearnedAdaptiveIndex/LAI.hpp"

#define DOWN_TIME_IN_SECONDS 5

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
    std::cout << "Learned Index Error Bound: " << error_bound << std::endl;
    std::cout << "Down Time: " << DOWN_TIME_IN_SECONDS << " seconds" <<std::endl;

    LAI <datatype_t> lai(dataset, inputDataSize, error_bound);
    double totalTimeElapsed = 0, queryTimeElapsed = 0, predictionTimeElapsed = 0, cumulativeQueryTimeOnly = 0;

    std::string query_folder = base_path + "queries_sel_random/";
    std::string real_query_file_path = query_folder + argv[3] + ".txt";
    std::ifstream query_file(real_query_file_path);
    std::string predicted_query_file_path = query_folder + "Predictions/" + argv[3] + ".txt";
    std::ifstream predictions_file(predicted_query_file_path);

    if (!query_file.is_open())
    {
        std::cout << "Query file is missing.. Abort.." <<endl;
        exit(0);
    }

    if (!predictions_file.is_open())
    {
        std::cout << "Prediction file is missing.. Abort.." <<endl;
        exit(0);
    }

    std::cout << "Query File Path = " << real_query_file_path << endl;
    std::cout << "Prediction File Path = " << predicted_query_file_path << endl;
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
            cumulativeQueryTimeOnly += time;
            query_num++;
            if(query_num == 1)
            {
	        std::cout << query_num << "," << cumulativeQueryTimeOnly/pow(10,9) << "," << cumulativeQueryTimeOnly/pow(10,9) << endl;
	        std::cout.flush();
            }
            
            //create index during down-time based on predictions
            if(query_num % subRuns == 0)
            {
	        int pred_query_count = 0;
                string line2;
                std::cout.flush();

                auto started = std::chrono::high_resolution_clock::now();        
                while(predictions_file >> line2)
                {
                    char *query_bounds = line2.data();
                    char *token = strtok(query_bounds, ",");
                    datatype_t query_num = atoi(token);
                    token = strtok(NULL, ",");

                    datatype_t lowVal = atoi(token);
                    if (lowVal >= inputDataSize)
                        lowVal = inputDataSize - 1;

                    token = strtok(NULL, ",");

                    datatype_t highVal = atoi(token);
                    if (highVal >= inputDataSize)
                        highVal = inputDataSize - 1;
                    
                    if (lowVal > highVal)
                        std::swap(lowVal, highVal);

                    if (lowVal == highVal)
                        continue;

                    std::pair<long,long> pos = lai.query(lowVal,highVal,false);
                    auto lowValPos = pos.first, highValPos = pos.second;
                    pred_query_count++;
                    if(pred_query_count == subRuns)
                        break;
                }
                auto done = std::chrono::high_resolution_clock::now();
                predictionTimeElapsed = (double)std::chrono::duration_cast<std::chrono::nanoseconds>(done-started).count();

                long timeToSleep = long(DOWN_TIME_IN_SECONDS * pow(10,9)) - predictionTimeElapsed;
                std::this_thread::sleep_for(std::chrono::nanoseconds(timeToSleep));

		done = std::chrono::high_resolution_clock::now();
                predictionTimeElapsed = (double)std::chrono::duration_cast<std::chrono::nanoseconds>(done-started).count();

                totalTimeElapsed += queryTimeElapsed + predictionTimeElapsed;
		std::cout << query_num << "," << cumulativeQueryTimeOnly/pow(10,9) << "," << totalTimeElapsed/pow(10,9) << endl;
		queryTimeElapsed = 0;
                     
            }
            
    }
    lai.getQueryStats();
}

