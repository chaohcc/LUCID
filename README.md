# LUCID
An Updatable and Concurrent Learned Index for Larger-than-Memory Data Management

## Getting started
LUCID can be used as a header-only library. To compile the program, you must enable at least the C++17 standard (e.g., use -std=c++17 in your compiler flags or set it in CMakeLists.txt). The minimum required CMake version is 3.21, as specified in the CMakeLists.txt file: cmake_minimum_required(VERSION 3.21)

To compile and run LUCID, use the [main.cpp] file included in the project.

There are some examples in [main.cpp]. e.g., -- test_merge_threshold: test LUCID on the blcok merge threshold, that when to merger a whole disk block. -- test_record: test the query workload with varying record sizes.  -- test_error: test LUCID under different error bound. -- test_blocksize: test the impact with different block sizes. 


## Dataset 

You can run this project on your own dataset by modifying the filepath in [main.cpp]. Your dataset should be formatted either as a binary file or as a text file (with one key per line). The datasets used in our paper are from [SOSD], [GRE], and [TLI]. Table 1 shows the details of datasets. 

#### <center>Tabel 1:  dataset description</center>
| dataset | description | duplicate | reference |
| :---:     |     :---:      |       :---:  |       :---:  |
books    | Amazon book sales popularity                                                                                    | No       | SOSD      |
fb        | Upsampled Facebook user ID                                                                                                                                 | No        |   SOSD     |
osm       | Uniformly sampled OpenStreetMap locations                                                                                                                  | No        |   SOSD    |
wiki      | Wikipedia article edit timestamps                                                                                                                          | Yes       |    SOSD    |
covid     | Uniformly sampled Tweet ID with tag COVID-19                                                                                                               | No        |    SNAM   |
genome    | Loci pairs in human chromosomes                  | No        |    Cell    |
stack     | Vote ID from Stackoverflow                                                                                                                                 | No       |     Stackoverflow   |
wise      | Partition key from the WISE data                                                                                                                           | No       |     AJ   |
libio     | Repository ID from libraries.io                                                                                                                            | No        |  Libraries.io.    |
history   | History node ID in OpenStreetMap                                                                                                                           | No        |    OpenStreetMap    |
planet    | Planet ID in OpenStreetMap                                                                                                                                 | No        |   OpenStreetMap    |
lognormal | Values generated according to a lognormal distribution, multiply $10^9$ and rounded down to the nearest integer | Yes       |    ALEX    | 

---
[SOSD] 1) Andreas Kipf and Ryan Marcus and Alexander van Renenr and others. SOSD: A Benchmark for Learned Indexes. MLSys@NeurIPS2019. 2) Ryan Marcus and Andreas Kipf and Alexander van Renen and others. Benchmarking Learned Indexes. VLDB 2020.



### References
[GRE] Chaichon Wongkham, Baotong Lu, Chris Liu, Zhicong Zhong, Eric Lo, and Tianzheng Wang. Are Updatable Learned Indexes Ready? . PVLDB,15(11): 3004 - 3017, 2022.

[GRE] Wongkham, Chaichon, Baotong Lu, Chris Liu, Zhicong Zhong, Eric Lo, and Tianzheng Wang. “Are Updatable Learned Indexes Ready?,” VLDB 2022.

[TLI] Sun, Zhaoyan, Xuanhe Zhou, and Guoliang Li. “Learned Index: A Comprehensive Experimental Evaluation.” VLDB 2023.

[SOSD] 1) Andreas Kipf and Ryan Marcus and Alexander van Renenr and others. SOSD: A Benchmark for Learned Indexes. MLSys@NeurIPS2019. 2) Ryan Marcus and Andreas Kipf and Alexander van Renen and others. Benchmarking Learned Indexes. VLDB 2020.

[SNAM] Christian E Lopez and Caleb Gallemore. 2021. An augmented multilingual Twitter dataset for studying the COVID-19 infodemic. Social Network Analysis and Mining 11, 1 (2021), 1–14.

[Cell] Suhas S.P. Rao and et al. 2014. A 3D Map of the Human Genome at Kilobase Resolution Reveals Principles of Chromatin Looping. Cell (2014).

[Stackoverflow] Stackoverflow. 2021. Vote ID. (2021). https://archive.org/download/stackexchange.

[Aj] Edward L Wright, Peter RM Eisenhardt, Amy K Mainzer, Michael E Ressler, Roc M Cutri, Thomas Jarrett, J Davy Kirkpatrick, Deborah Padgett, Robert S McMillan, Michael Skrutskie, et al. 2010. The Wide-field Infrared Survey Explore (WISE): mission description and initial on-orbit performance. The Astronomical Journal 140, 6 (2010), 1868.

[Libraries.io.] Libraries.io. 2017. Repository ID. (2017). https://libraries.io/data. 

[OpenStreetMap] Google Cloud. 2017. OpenStreetMap. (2017). https://console.cloud.google.com/marketplace/details/openstreetmap/geo-openstreetmap

[ALEX] Jialin Ding, Umar Farooq Minhas, Jia Yu, Chi Wang, Jaeyoung Do, Yinan Li,Hantian Zhang, Badrish Chandramouli, Johannes Gehrke, Donald Kossmann, et al. 2020. ALEX: an updatable adaptive learned index. Proceedings of the 2020 ACM SIGMOD International Conference on Management of Data, 969–984.



