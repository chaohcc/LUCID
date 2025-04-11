//
// Created by chaochao on 2023/2/14.
//

#ifndef TOTALY_REBUILD_FILM_WORKLOAD_H
#define TOTALY_REBUILD_FILM_WORKLOAD_H


#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <algorithm>
#include <atomic>
#include <memory>
#include <random>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <thread>
#include <ctime>
#include <algorithm>

//#include "parallel_sort.h"

#include "flags.h"
#include "utils.h"

template<typename KEY_TYPE, typename PAYLOAD_TYPE>
class Workload{
public:
    enum Operation {
        READ = 0, INSERT, DELETE, SCAN, UPDATE, DIFF
    };

    double read_ratio = 1;
    double insert_ratio = 0;
    double delete_ratio = 0;
    double update_ratio = 0;
    double scan_ratio = 0;
    size_t scan_num = 100;
    size_t operations_num;
    long long table_size = -1;
    size_t init_table_size;
    double init_table_ratio;
    double del_table_ratio;
    double buffer_ratio;
    size_t thread_num = 1;
    int batch_num = 0;
    size_t delete_counter = 0;
    std::vector <std::string> all_index_type;
    std::vector <std::string> all_thread_num;
//    std::string index_type;
    std::string keys_file_type;
    std::string sample_distribution;
    bool latency_sample = false;
    double latency_sample_ratio = 0.01;
    int error_bound;
    std::string output_path;
    size_t random_seed;
    double zipf_factor = 0.75;
    bool memory_record;
    bool dataset_statistic;
    bool data_shift = false;

    std::vector <KEY_TYPE> init_keys;
    KEY_TYPE *keys;
    std::pair <KEY_TYPE, PAYLOAD_TYPE> *init_key_values;
//    std::vector <std::pair<Operation, KEY_TYPE>> operations;
    std::vector <KEY_TYPE> read_ops;
    std::vector <KEY_TYPE> update_ops;
    std::vector <KEY_TYPE> insert_ops;
    std::vector <KEY_TYPE> scan_ops;   // the start key of scan range (with scan length fixed, or with various length)
    std::vector <KEY_TYPE> delete_ops;
    KEY_TYPE** range_ops = nullptr;    // the end key of scan range (the scan range with various length)
    std::mt19937 gen;

    Workload() {
    }


    inline void parse_args(int argc, char **argv) {
        auto flags = parse_flags(argc, argv);
        keys_file_path = get_required(flags, "keys_file"); // required
        keys_file_type = get_with_default(flags, "keys_file_type", "binary");
        read_ratio = stod(get_required(flags, "read")); // required
        insert_ratio = stod(get_with_default(flags, "insert", "0")); // required
        delete_ratio = stod(get_with_default(flags, "delete", "0"));
        update_ratio = stod(get_with_default(flags, "update", "0"));
        buffer_ratio = stod(get_with_default(flags, "buffer", "0.3"));
        scan_ratio = stod(get_with_default(flags, "scan", "0"));
        scan_num = stoi(get_with_default(flags, "scan_num", "100"));
        operations_num = stoi(get_with_default(flags, "operations_num", "1000000000")); // required
        table_size = stoi(get_with_default(flags, "table_size", "-1"));
        init_table_ratio = stod(get_with_default(flags, "init_table_ratio", "1"));
        del_table_ratio = stod(get_with_default(flags, "del_table_ratio", "0.5"));
        init_table_size = 0;
        all_thread_num = get_comma_separated(flags, "thread_num"); // required
        all_index_type = get_comma_separated(flags, "index"); // required
        sample_distribution = get_with_default(flags, "sample_distribution", "uniform");
        latency_sample = get_boolean_flag(flags, "latency_sample");
        latency_sample_ratio = stod(get_with_default(flags, "latency_sample_ratio", "0.01"));
        error_bound = stoi(get_with_default(flags, "error_bound", "64"));
        output_path = get_with_default(flags, "output_path", "./out.csv");
        random_seed = std::random_device{}();
        random_seed = stoul(get_with_default(flags, "seed", "412773089"));
//        random_seed = stoul(get_with_default(flags, "seed", std::to_string(random_seed)));
        gen.seed(random_seed);
//        gen.seed(0);
        memory_record = get_boolean_flag(flags, "memory");
        dataset_statistic = get_boolean_flag(flags, "dataset_statistic");
        data_shift = get_boolean_flag(flags, "data_shift");
        batch_num = stoi(get_with_default(flags,"batch_num","11"));

        COUT_THIS("[micro] Read:Insert:Update:Scan:Delete= " << read_ratio << ":" << insert_ratio << ":" << update_ratio << ":"
                                                             << scan_ratio << ":" << delete_ratio);
        double ratio_sum = read_ratio + insert_ratio + delete_ratio + update_ratio + scan_ratio;
        double insert_delete = insert_ratio + delete_ratio;
//        INVARIANT(insert_delete == insert_ratio || insert_delete == delete_ratio);   //这说明在一个workload 中不能同时存在insert和delete  ； 注释掉以后 在这个workload 中可以同时存在insert和delete
        INVARIANT(ratio_sum > 0.9999 && ratio_sum < 1.0001);  // avoid precision lost
        INVARIANT(sample_distribution == "zipf" || sample_distribution == "uniform" ||
                  sample_distribution == "random" || sample_distribution == "zipfrandom" || sample_distribution == "hotspot");

    }

    void generate_operations(KEY_TYPE *keys) {
        size_t range_num = 0;
        // prepare operations
        std::cout << "i need You, my Lord, start generate operations, sample keys according to " << sample_distribution << std::endl;
        read_ops.reserve(operations_num * read_ratio +200);
        insert_ops.reserve(operations_num * insert_ratio +200);
        update_ops.reserve(operations_num * update_ratio +200);
        scan_ops.reserve(operations_num * scan_ratio +200);
        delete_ops.reserve(operations_num * delete_ratio +200);
        KEY_TYPE *sample_ptr = nullptr;
        if (sample_distribution == "uniform") {
            sample_ptr = get_search_keys(&init_keys[0], init_table_size, operations_num, &random_seed);
        } else if (sample_distribution == "zipf") {
            sample_ptr = get_search_keys_zipf(&init_keys[0], init_table_size, operations_num, zipf_factor);
        }

        // generate operations(read, insert, update, scan)
        COUT_THIS("generate operations.");
        std::uniform_real_distribution<> ratio_dis(0, 1);
        size_t sample_counter = 0, insert_counter = init_table_size;
        size_t delete_counter = table_size * (1 - del_table_ratio);

        if (data_shift) {
            size_t rest_key_num = table_size - init_table_size;
            if(rest_key_num > 0) {
                std::sort(keys + init_table_size, keys + table_size);
//                tbb::parallel_sort(keys + init_table_size, keys + table_size);
                std::random_shuffle(keys + init_table_size, keys + table_size);
            }
        }

        size_t temp_counter = 0;
        for (size_t i = 0; i < operations_num; ++i) {
            auto prob = ratio_dis(gen);
            if (prob < read_ratio) {
                // if (temp_counter >= table_size) {
                //     operations_num = i;
                //     break;
                // }
                // operations.push_back(std::pair<Operation, KEY_TYPE>(READ, keys[temp_counter++]));
                read_ops.push_back(sample_ptr[sample_counter++]);
            } else if (prob < read_ratio + insert_ratio) {
                if (insert_counter >= table_size) {
                    operations_num = i;
                    break;
                }
                insert_ops.push_back( keys[insert_counter++]);
            } else if (prob < read_ratio + insert_ratio + update_ratio) {
                update_ops.push_back(sample_ptr[sample_counter++]);
            } else if (prob < read_ratio + insert_ratio + update_ratio + scan_ratio) {
                scan_ops.push_back(sample_ptr[sample_counter++]);
            } else {
                if (delete_counter >= table_size) {
                    operations_num = i;
                    break;
                }
                delete_ops.push_back( keys[delete_counter++]);
                // operations.push_back(std::pair<Operation, KEY_TYPE>(DELETE, sample_ptr[sample_counter++]));
            }
        }
        long batch_operation_num = insert_ops.size()+read_ops.size()+update_ops.size() + scan_ops.size()+delete_ops.size();
        COUT_VAR(batch_operation_num);

        delete[] sample_ptr;
    }

    double generate_operations(KEY_TYPE *keys, unsigned long pos_shift) {
        // prepare operations
        struct timeval start_time,end_time;
        gettimeofday(&start_time, NULL);
        size_t range_num = 0;
        std::cout << "i need You, my Lord, start generate operations, sample keys according to " << sample_distribution << std::endl;
//        operations.reserve(operations_num);
        read_ops.reserve(operations_num * read_ratio +200);
        insert_ops.reserve(operations_num * insert_ratio +200);
        update_ops.reserve(operations_num * update_ratio +200);
        scan_ops.reserve(operations_num * scan_ratio +200);
        delete_ops.reserve(operations_num * delete_ratio +200);
//        long batch_operation_num1 = insert_ops.size()+read_ops.size()+update_ops.size() + scan_ops.size()+delete_ops.size();
//        COUT_VAR(batch_operation_num1);
        KEY_TYPE *sample_ptr = nullptr;

        if (pos_shift > init_table_size) {
            // renew init_keys, those keys have been inserted into index
            if (read_ratio > 0 || scan_ratio > 0){
                init_keys.insert(init_keys.end(), keys + init_table_size, keys + pos_shift);
//            tbb::parallel_sort(init_keys.begin(), init_keys.end());
                std::sort(init_keys.begin(), init_keys.end());
            }

            init_table_size = pos_shift;
        }

        if (sample_distribution == "uniform") {
            sample_ptr = get_search_keys(&init_keys[0], init_table_size, operations_num, &random_seed);
        } else if (sample_distribution == "zipf") {
            sample_ptr = get_search_keys_zipf(&init_keys[0], init_table_size, operations_num, zipf_factor);
        }else if (sample_distribution == "random") {
            sample_ptr = get_search_keys(&init_keys[0], init_table_size, operations_num);
        } else if (sample_distribution == "zipfrandom") {
            sample_ptr = get_search_keys_scrambledzipf(&init_keys[0], init_table_size, operations_num);
        } else if (sample_distribution == "hotspot") {
            sample_ptr = get_search_keys_hotspot(&init_keys[0], init_table_size, operations_num, zipf_factor);
        } else {
            std::cerr << "--lookup_distribution must be either 'hotspot', 'randomzipf','random' or 'zipf'"
                      << std::endl;
            exit;
        }
        size_t random_seed = 521959958;
        auto random_sample_ptr = get_delete_keys(&init_keys[0], init_table_size, operations_num*9, &random_seed);

        if (data_shift) {
            size_t rest_key_num = table_size - init_table_size;
            if(rest_key_num > 0) {
                std::sort(keys + init_table_size, keys + table_size);
//                tbb::parallel_sort(keys + init_table_size, keys + table_size);
                std::random_shuffle(keys + init_table_size, keys + table_size);
            }
        }


        // generate operations(read, insert, update, scan)
        COUT_THIS("generate operations.");
        std::uniform_real_distribution<> ratio_dis(0, 1);
        size_t sample_counter = 0, insert_counter = init_table_size;
//        size_t delete_counter = init_table_size * (1 - del_table_ratio);
//        size_t delete_counter = std::ceil(((batch_no-1)*delete_ratio)*operations_num) + 2;
        for (size_t i = 0; i < operations_num; ++i) {
            auto prob = ratio_dis(gen);
            if (prob < read_ratio) {
                // if (temp_counter >= table_size) {
                //     operations_num = i;
                //     break;
                // }
                // operations.push_back(std::pair<Operation, KEY_TYPE>(READ, keys[temp_counter++]));
//               if (sample_ptr[sample_counter] == 62129282194 ){
//                  std::cout<< "Jesus, please come!, 62129282194, sample_counter is " << sample_counter<<std::endl;
//               }
                read_ops.push_back(sample_ptr[sample_counter++]);
            } else if (prob < read_ratio + insert_ratio) {
                if (insert_counter >= table_size) {
                    operations_num = i;
                    break;
                }
//                if (keys[insert_counter] == 62129282194 ){
//                  std::cout<< "Jesus, please come!, 62129282194, insert_counter is " << insert_counter<<std::endl;
//                }
                insert_ops.push_back( keys[insert_counter++]);
            } else if (prob < read_ratio + insert_ratio + update_ratio) {
                update_ops.push_back(sample_ptr[sample_counter++]);
            } else if (prob < read_ratio + insert_ratio + update_ratio + scan_ratio) {
                scan_ops.push_back(sample_ptr[sample_counter++]);  // scan with a fixed length
                range_num++;
            } else {
                if (delete_counter >= table_size) {
                    operations_num = i;
                    break;
                }
                delete_ops.push_back( random_sample_ptr[delete_counter++]);
                // operations.push_back(std::pair<Operation, KEY_TYPE>(DELETE, sample_ptr[sample_counter++]));
            }
        }

        long batch_operation_num = insert_ops.size()+read_ops.size()+update_ops.size() + scan_ops.size()+delete_ops.size();
        COUT_VAR(batch_operation_num);
        init_table_size = pos_shift;
        if (range_num){
            // generate range query according to the workload distribution
            if (sample_distribution == "random") {
                range_ops = get_search_ranges(init_keys, init_table_size, range_num);
            } else if (sample_distribution == "zipfrandom") {
                range_ops = get_search_ranges_scrambledzipf(init_keys, init_table_size, range_num);
            } else if (sample_distribution == "zipf") {
                range_ops = get_search_ranges_zipf(init_keys, init_table_size, range_num, zipf_factor,0,scan_num);
            } else if (sample_distribution == "hotspot") {
                range_ops = get_search_ranges_hotspot(init_keys, init_table_size,range_num, zipf_factor);
            } else {
                std::cerr << "--lookup_distribution must be either 'hotspot', 'zipfrandom','random' or 'zipf'"
                          << std::endl;
            }
        }

        delete[] sample_ptr;

        gettimeofday(&end_time, NULL);
        double workload_elapsed_time = (end_time.tv_sec - start_time.tv_sec) +
                                       (double) (end_time.tv_usec - start_time.tv_usec) / 1000000.0;
        return workload_elapsed_time;
    }

    KEY_TYPE *load_keys(long long test_size) {
        // Read keys from file
        COUT_THIS("Reading data from file.");

        if (table_size > 0) keys = new KEY_TYPE[table_size];
        if (keys_file_type == "binary") {
            table_size = load_binary_data(keys, table_size, keys_file_path);
            if (table_size <= 0) {
                COUT_THIS("Could not open key file, please check the path of key file.");
                exit(0);
            }
        } else if (keys_file_type == "text") {
            table_size = load_text_data(keys, table_size, keys_file_path);
            if (table_size <= 0) {
                COUT_THIS("Could not open key file, please check the path of key file.");
                exit(0);
            }
        } else {
            COUT_THIS("Could not open key file, please check the path of key file.");
            exit(0);
        }

        if (!data_shift) {
//            tbb::parallel_sort(keys, keys + table_size);
            if (keys_file_type != "text")
                std::sort(keys, keys + table_size);
            auto last = std::unique(keys, keys + table_size);
            table_size = last - keys;
            if (table_size > test_size){

                KEY_TYPE* nn = (KEY_TYPE*)realloc(keys, test_size*sizeof(KEY_TYPE));// 截断数组
                table_size = test_size;
            }
            std::shuffle(keys, keys + table_size, gen);
        }

//        init_table_size = init_table_ratio * table_size;
        std::cout << "Table size is " << table_size << ", Initial table size is " << init_table_size << std::endl;

        for (auto j = 0; j < 10; j++) {
            std::cout << keys[j] << " ";
        }
        std::cout << std::endl;

        // prepare data
        COUT_THIS("prepare init keys.");
        init_keys.resize(init_table_size);
        for (size_t i = 0; i < init_table_size; ++i) {
            init_keys[i] = (keys[i]);
        }
//        tbb::parallel_sort(init_keys.begin(), init_keys.end());
        std::sort(init_keys.begin(), init_keys.end());
//        auto max_key = init_keys.back(); //161142159158

//        init_key_values = new std::pair<KEY_TYPE, PAYLOAD_TYPE>[init_keys.size()];
//#pragma omp parallel for num_threads(thread_num)
//        for (int i = 0; i < init_keys.size(); i++) {
//            init_key_values[i].first = init_keys[i];
////            init_key_values[i].second = 123456789;
//        }
        COUT_VAR(table_size);
//        COUT_VAR(init_keys.size());

        return keys;
    }

    std::string keys_file_path;
};



#endif //TOTALY_REBUILD_FILM_WORKLOAD_H
