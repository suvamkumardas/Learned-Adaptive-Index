# Learned-Adaptive-Index

There are two .cpp files:
- test_lai_with_predictions.cpp (For LAI with predictions enabled)
- test_lai.cpp (For LAI with predictions disabled)

The data files must be stored at "data_and_queries/100_million/". As of now the data file is absent due to space restrictions.

The query workloads are stored in the directory "data_and_queries/100_million/queries_sel_random/". Each workload has a separate ".txt" file for storing the queries. For the sake of experiments, the predicted query files are stored at "data_and_queries/100_million/queries_sel_random/Predictions/" directory. For each query workload, there is a predicted workload file present.

1. To run LAI with predictions enabled, compile test_lai_with_predictions.cpp as follows:
g++ -O3 test_lai_with_predictions.cpp -o test_lai_with_predictions

2. To run LAI with predictions disabled, compile test_lai.cpp as follows:
g++ -O3 test_lai.cpp -o test_lai

3. To run LAI execute:
<exe-file> <input-data-file> <data-size> <workload> <subRuns-to-report-cumulative-execution-times>
e.g.
	- ./test_lai_with_predictions data.txt 100000000 Random 1000 (To run LAI with predictions on Random workload)
	- ./test_lai data.txt 100000000 SeqAlt 1000 (To run LAI without predictions on SeqAlt workload)
	
4. Please note that, to execute LAI the data must be read from a file stored in "data_and_queries/100_million/" directory. To run other workloads, the name of the workload file should be provided as a command line argument, as shown above, without the ".txt" extension.

