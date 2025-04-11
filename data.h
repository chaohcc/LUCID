//
// Created by chaochao on 2021/12/17.
//

#ifndef DATA_H
#define DATA_H

#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <fstream>
#include<ctime>
#include<iostream>
#include <string>
#include<fstream>
#include <vector>
#include <random>
#include "zipf.h"
#include "utils.h"


template<class T>
long long load_binary_data(T *&data, long long length, const std::string &file_path) {
    // open key file
    std::ifstream is(file_path.c_str(), std::ios::binary | std::ios::in);
    if (!is.is_open()) {
        return 0;
    }

    std::cout << "The path of data file: " << file_path << std::endl;

    // read the number of keys
    T max_size;
    is.read(reinterpret_cast<char*>(&max_size), sizeof(T));


    std::cout << "the max size of data file:" << max_size << std::endl;

    // create array
//    if(length < 0 || length > max_size) length = max_size;   // comented by chao
    if(length < 0  ) length = max_size;   // chao changed
    data = new T[length];

    // read keys
    is.read(reinterpret_cast<char *>(data), std::streamsize(length * sizeof(T)));
    is.close();
    return length;
}


template <class T>
long load_text_data(T *&array, int length, const std::string& file_path) {
    std::ifstream is(file_path.c_str());
    if (!is.is_open()) {
        return false;
    }
    long i = 0;
    std::string str;
    std::vector<T> temp_keys;
    temp_keys.reserve(2000000000);
    double a;
    auto aj = sizeof(T);
    while (std::getline(is, str) && (i < length )) {
        std::istringstream ss(str);
        ss >> a;
//        array[i] = a;
        temp_keys.push_back(a);
//        ss >> array[i];
        i++;
    }
    array = new T[temp_keys.size()];
    if (!temp_keys.empty()){
        memcpy(array, &temp_keys[0],temp_keys.size()*sizeof(T));
    }
    is.close();
    std::vector<T>().swap(temp_keys);
    return i;

}



/// load the 'dset' dataset, the number of loaded keys is 'numkeys'
template<typename type_key>
bool loaddata(type_key data[], std::string dset, size_t numkeys) {
    std::string filepath;
    //std::vector<type_key> keys (numkeys);

    if (!dset.compare("books")) {
        filepath = "/home/wamdm/chaohong/clionDir/FeasFearCPP/dataset/books.txt";
        load_text_data(data, numkeys, filepath);
    } else if (!dset.compare("random0.5")) {
        filepath = "/home/wamdm/chaohong/clionDir/FeasFearCPP/dataset/random0.5_83886080.txt";
        load_text_data(data, numkeys, filepath);
    } else if (!dset.compare("wiki_ts")) {
        filepath = "/home/wamdm/chaohong/clionDir/FeasFearCPP/dataset/wiki_ts_100000000.txt";
        load_text_data(data, numkeys, filepath);
    } else if (!dset.compare("astro_ra")) {
        filepath = "/home/wamdm/chaohong/clionDir/FeasFearCPP/dataset/astro_ra_18_21_89523861.txt";
        load_text_data(data, numkeys, filepath);
    }
    else if (!dset.compare("lognormal")){
        filepath = "/home/wamdm/chaohong/dataAchive/ALEX_1D/lognormal-190M.bin.data";
        load_binary_data(data, numkeys, filepath);
    }
    else if (!dset.compare("longlat")){
        filepath = "/home/wamdm/chaohong/dataAchive/ALEX_1D/longlat-200M.bin.data";
        load_binary_data(data, numkeys, filepath);
    }
    else if (!dset.compare("longitudes")){
        filepath = "/home/wamdm/chaohong/dataAchive/ALEX_1D/longitudes-200M.bin.data";
        load_binary_data(data, numkeys, filepath);
    }
    else if (!dset.compare("ycsb")){
        filepath = "/home/wamdm/chaohong/dataAchive/ALEX_1D/ycsb-200M.bin.data";
        load_binary_data(data, numkeys, filepath);
    }
    else {
        printf("data name is wrong");
        exit(0);
    }
    return true;
}

template<typename type_key>
bool loaddata(type_key data[], const std::string &keys_file_path, std::string keys_file_type, size_t numkeys) {
    //std::vector<type_key> keys (numkeys);
    if (keys_file_type == "binary"){
        load_binary_data(data, numkeys, keys_file_path);
    }
    else if (keys_file_type == "text") {
        load_text_data(data, numkeys, keys_file_path);
    }else{
        COUT_THIS("Could not open key file, please check the path of key file.");
        exit(0);
    }
}

std::vector<std::string> split(std::string str, std::string pattern)
{
    std::string::size_type pos;
    std::vector<std::string> result;
    str += pattern;//扩展字符串以方便操作
    int size = str.size();
    for (int i = 0; i < size; i++)
    {
        pos = str.find(pattern, i);
        if (pos < size)
        {
            std::string s = str.substr(i, pos - i);
            result.push_back(s);
            i = pos + pattern.size() - 1;
        }
    }
    return result;
}


template <class T>
T* get_search_keys(T* array, int num_keys, int num_searches, size_t *seed = nullptr) {
//    std::mt19937_64 gen(std::random_device{}());
    std::mt19937_64 gen(9958);
    std::uniform_int_distribution<int> dis(0, num_keys - 1);
    auto* keys = new T[num_searches];
    for (int i = 0; i < num_searches; i++) {
        int pos = dis(gen);
        keys[i] = array[pos];
    }
    return keys;
}

template <class T>
T* get_delete_keys(T* array, int num_keys, int num_searches, size_t *seed = nullptr) {
//    std::mt19937_64 gen(std::random_device{}());
    if (num_searches > num_keys) {
        throw std::invalid_argument("num_searches cannot be greater than num_keys");
    }

    std::mt19937_64 gen(seed ? *seed : 9958);
    T* keys = new T[num_searches];
//    创建索引数组
    std::vector<int> indices(num_keys);
//    串行初始化
//    for (int i = 0; i < num_keys; i++) indices[i] = i;
    // 打乱 (全部shuffle)  速度较慢
//    std::shuffle(indices.begin(), indices.end(), gen);


    // 串行：Partial Fisher-Yates Shuffle: 仅随机选择前 num_searches 个
//    for (int i = 0; i < num_searches; i++) {
//        std::uniform_int_distribution<int> dis(i, num_keys - 1);
//        int swap_idx = dis(gen);
//        std::swap(indices[i], indices[swap_idx]);
//    }
//    // 串行 选择前 num_searches 个元素
//    for (int i = 0; i < num_searches; i++) {
//        keys[i] = array[indices[i]];
//    }

    // **并行初始化索引数组**
#pragma omp parallel for
    for (int i = 0; i < num_keys; i++) {
        indices[i] = i;
    }

    // **并行部分 Fisher-Yates 洗牌**
#pragma omp parallel for
    for (int i = 0; i < num_searches; i++) {
        std::uniform_int_distribution<int> dis(i, num_keys - 1);
        int swap_idx = dis(gen);
        std::swap(indices[i], indices[swap_idx]);
    }
//
#pragma omp parallel for
    for (int i = 0; i < num_searches; i++) {
        keys[i] = array[indices[i]];
    }

    return keys;
}

template <class T>
T* get_search_keys_scrambledzipf(T* array, unsigned long int num_keys, unsigned int num_searches, size_t *seed = nullptr) {
    auto* keys = new T[num_searches];
    ScrambledZipfianGenerator zipf_gen(num_keys);
    for (unsigned int i = 0; i < num_searches; i++) {
        unsigned int pos = zipf_gen.nextValue();
        keys[i] = array[pos];
    }
    return keys;
}

template <class T>
T* get_search_keys_zipf(T* array, unsigned long int num_keys, unsigned long int num_searches, double zipf_factor) {

    auto* keys = new T[num_searches];
    zipf_distribution<long int, double > zipf_gen(num_keys,zipf_factor);
//    std::mt19937 gen_;
    std::mt19937 gen_(0);
    std::vector<T> zipf_res;
    for (unsigned int i = 0; i< num_searches; i ++){
        unsigned int pos = num_keys - zipf_gen.operator()(gen_);
        keys[i] = array[pos];
    }
    return keys;
}

/*
template <class T>
T* get_search_keys_hotspot(std::vector<T> array, unsigned long int num_keys, unsigned long int num_searches, double hotratio=0.2,double accessratio = 0.9) {
    unsigned int hotspotlen = num_keys * hotratio;
    unsigned int hotqueryn = num_searches*accessratio;
    unsigned int randomqueryn = num_searches - hotqueryn;
    auto* keys = new T[num_searches];
    zipf_distribution<long int, double > zipf_gen(num_keys,0.75);
    std::mt19937_64 gen_random(std::random_device{}());
    std::mt19937 gen_;
    std::vector<T> hospot_res;
    unsigned int hot_start = num_keys - zipf_gen.operator()(gen_);
    while (num_keys-hot_start < hotspotlen)
        hot_start = hot_start-hotspotlen;
    std::uniform_int_distribution<int> dis1(hot_start, hot_start+hotspotlen);
    for (unsigned int i = 0; i< hotqueryn; i ++){
        unsigned int pos  = dis1(gen_random);
        keys[i] = array[pos];
    }

    std::uniform_int_distribution<int> dis2(0, num_keys - 1);
    for (unsigned int i = 0;i<randomqueryn;i++){
        unsigned int pos = dis2(gen_random);
        keys[hotqueryn+i] = array[pos];
    }
    return keys;
}
*/

template <class T>
T* get_search_keys_hotspot(T* array, unsigned long int num_keys, unsigned long int num_searches, double hotratio=0.2,double accessratio = 0.9) {
    unsigned int hotspotlen = num_keys * hotratio;
    unsigned int hotqueryn = num_searches*accessratio;
    unsigned int randomqueryn = num_searches - hotqueryn;
    auto* keys = new T[num_searches];
    zipf_distribution<long int, double > zipf_gen(num_keys,0.75);
    std::mt19937_64 gen_random(std::random_device{}());
//    std::mt19937_64 gen_random(0);
    std::mt19937 gen_;
//    std::mt19937 gen_(0);
    std::vector<T> hospot_res;
    unsigned int hot_start = num_keys - zipf_gen.operator()(gen_);
    while (num_keys-hot_start < hotspotlen)
        hot_start = hot_start-hotspotlen;
    std::uniform_int_distribution<int> dis1(hot_start, hot_start+hotspotlen);
    for (unsigned int i = 0; i< hotqueryn; i ++){
        unsigned int pos  = dis1(gen_random);
        keys[i] = array[pos];
    }

    std::uniform_int_distribution<int> dis2(0, num_keys - 1);
    for (unsigned int i = 0;i<randomqueryn;i++){
        unsigned int pos = dis2(gen_random);
        keys[hotqueryn+i] = array[pos];
    }
    return keys;
}


// generate the range queries
template <class T>
T** get_search_ranges(std::vector<T> array, int num_keys, int num_searches,int minlen = 1,int maxlen =100) {
//    std::mt19937_64 gen(std::random_device{}());
    std::mt19937_64 gen(9958);
    std::uniform_int_distribution<int> dis(0, num_keys - maxlen);
//    std::mt19937_64 gen_random(std::random_device{}());
    std::mt19937_64 gen_random(1314);
    std::uniform_int_distribution<int> disrange(minlen, maxlen);
    auto* ranges = new T*[num_searches];
    for (int i = 0; i < num_searches; i++) {
        unsigned int pos = dis(gen);
        ranges[i] = new T[2];
        ranges[i][0] = array[pos];
        ranges[i][1] = array[pos+disrange(gen_random)];
    }
    return ranges;
}


template <class T>
T** get_search_ranges_zipf(std::vector<T> array, unsigned long int num_keys, unsigned long int num_searches, double zipf_factor,unsigned int minlen = 0,unsigned int maxlen=100) {

    auto* ranges = new T*[num_searches];
    long int upper = num_keys-maxlen;
    zipf_distribution<long int, double > zipf_gen(upper,zipf_factor);
//    std::mt19937 gen_;
    std::mt19937 gen_(0);
    std::uniform_int_distribution<int> disrange(minlen, maxlen);
//    std::mt19937_64 gen_random(std::random_device{}());
    std::mt19937_64 gen_random(0);
    for (int i = 0; i< num_searches; i ++){
        unsigned int pos = upper - zipf_gen.operator()(gen_);
        ranges[i] = new T[2];
        ranges[i][0] = array[pos];
        unsigned int pos2 = pos+disrange(gen_random);
//        if (pos2>num_keys-1){
//            std::cout<< "my Lord, i need You!"<< std::endl;
//        }
        ranges[i][1] = array[pos2];
    }
    return ranges;
}

template <class T>
T** get_search_ranges_scrambledzipf(std::vector<T> array, unsigned long int num_keys, unsigned int num_searches,int minlen = 0,int maxlen=100) {

    auto* ranges = new T*[num_searches];
    ScrambledZipfianGenerator zipf_gen(num_keys-maxlen);
    std::mt19937 gen_;
    std::uniform_int_distribution<int> disrange(minlen, maxlen);
    std::mt19937_64 gen_random(std::random_device{}());
    for (unsigned int i = 0; i< num_searches; i ++){
        unsigned int pos = zipf_gen.nextValue();
        ranges[i] = new T[2];
        ranges[i][0] = array[pos];
        ranges[i][1] = array[pos+disrange(gen_random)];
    }
    return ranges;
}


template <class T>
T** get_search_ranges_hotspot(std::vector<T> array, unsigned long int num_keys, unsigned long int num_searches, double hotratio=0.2,double accessratio = 0.9,int minlen = 0,int maxlen=100) {
    unsigned int hotspotlen = num_keys * hotratio;
    unsigned int hotqueryn = num_searches*accessratio;
    unsigned int randomqueryn = num_searches - hotqueryn;
    auto* ranges = new T*[num_searches];
    zipf_distribution<long int, double > zipf_gen(num_keys,0.5);
    std::mt19937_64 gen_random(std::random_device{}());
    std::mt19937 gen_;

    unsigned int hot_start = num_keys - zipf_gen.operator()(gen_);
    while (num_keys-hot_start< hotspotlen )
        hot_start = hot_start-hotspotlen;
    while (hot_start+hotspotlen+maxlen> num_keys )
        hot_start = hot_start-2*maxlen;
    std::uniform_int_distribution<int> dis1(hot_start, hot_start+hotspotlen);
    std::uniform_int_distribution<int> disrange(minlen, maxlen);
    for (unsigned int i = 0; i< hotqueryn; i ++){
        unsigned int pos  = dis1(gen_random);
        if (pos<hot_start )
            std::cout<<"i need You, my lovely Lord"<<std::endl;
        if (pos>hot_start+hotspotlen)
            std::cout<<"i need You, my lovely Lord"<<std::endl;
        ranges[i] = new T[2];
        ranges[i][0] = array[pos];
        ranges[i][1] = array[pos+disrange(gen_random)];

    }

    std::uniform_int_distribution<int> dis2(0, num_keys - maxlen);
    for (unsigned int i = 0;i<randomqueryn;i++){
        unsigned int pos = dis2(gen_random);
        ranges[hotqueryn+i] = new T[2];
        ranges[hotqueryn+i][0] = array[pos];
        ranges[hotqueryn+i][1] = array[pos+disrange(gen_random)];
    }
    return ranges;
}




#endif //DATA_H
