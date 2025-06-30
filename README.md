# LUCID
An Updatable and Concurrent Learned Index for Larger-than-Memory Data Management

## Getting started
LUCID can be used as a header-only library. To compile the program, you must enable at least the C++17 standard (e.g., use -std=c++17 in your compiler flags or set it in CMakeLists.txt). The minimum required CMake version is 3.21, as specified in the CMakeLists.txt file: cmake_minimum_required(VERSION 3.21)

To compile and run LUCID, use the [main.cpp] file included in the project.

There are some examples in [main.cpp]. e.g., -- test_merge_threshold: test LUCID on the blcok merge threshold, that when to merger a whole disk block. -- test_record: test the query workload with varying record sizes.  -- test_error: test LUCID under different error bound. -- test_blocksize: test the impact with different block sizes. 


## Dataset 

You can run this project on your own dataset by modifying the filepath in [main.cpp]. Your dataset should be formatted either as a binary file or as a text file (with one key per line). The datasets used in our paper are from [SOSD], [GRE], and [TLI].

[SOSD] 1) Andreas Kipf and Ryan Marcus and Alexander van Renenr and others. SOSD: A Benchmark for Learned Indexes. MLSys@NeurIPS2019. 2) Ryan Marcus and Andreas Kipf and Alexander van Renen and others. Benchmarking Learned Indexes. VLDB 2020
[GRE] Wongkham, Chaichon, Baotong Lu, Chris Liu, Zhicong Zhong, Eric Lo, and Tianzheng Wang. “Are Updatable Learned Indexes Ready?,” VLDB 2022.
[TLI] Sun, Zhaoyan, Xuanhe Zhou, and Guoliang Li. “Learned Index: A Comprehensive Experimental Evaluation.” VLDB 2023.



