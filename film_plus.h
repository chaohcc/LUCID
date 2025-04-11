//
// Created by CCMa on 2022/12/5.
//

#ifndef PLUSFILM_FILM_PLUS_H
#define PLUSFILM_FILM_PLUS_H

#include <stdlib.h>
#include <queue>
#include <sys/types.h>


#include <algorithm>
#include <array>
#include <chrono>
#include <cstring>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <set>
#include <string>
#include <utility>
#include <vector>
#ifdef _MSC_VER
#include <intrin.h>
#else
#include <x86intrin.h>
#endif
#include <bitset>
#include <cassert>
#ifdef _WIN32
#include <intrin.h>
#include <limits.h>
typedef unsigned __int32 uint32_t;
#else
#include <stdint.h>
#endif

#ifdef _MSC_VER
#define forceinline __forceinline
#elif defined(__GNUC__)
#define forceinline inline __attribute__((__always_inline__))
#elif defined(__CLANG__)
#if __has_attribute(__always_inline__)
#define forceinline inline __attribute__((__always_inline__))
#else
#define forceinline inline
#endif
#else
#define forceinline inline
#endif

// Whether we store key and payload arrays separately in data nodes
// By default, we store them separately
#define ALEX_DATA_NODE_SEP_ARRAYS 1

#if ALEX_DATA_NODE_SEP_ARRAYS
#define ALEX_DATA_NODE_KEY_AT(i) key_slots_[i]
//#define ALEX_DATA_NODE_PAYLOAD_AT(i) payload_slots_[i]
#define ALEX_DATA_NODE_PAYLOAD_AT(i) pointer_slots_[i].second
#define ALEX_DATA_NODE_FLAG_AT(i) pointer_slots_[i].first
#else
#define ALEX_DATA_NODE_KEY_AT(i) data_slots_[i].first
#define ALEX_DATA_NODE_PAYLOAD_AT(i) data_slots_[i].second
#endif

// *** Required Headers

#include <algorithm>
#include <functional>
#include <iostream>
#include <fstream>
#include <memory>
#include <cstddef>
#include <cassert>
#include <map>
#include <string>
#include <vector>
#include <cstdio>
#include <sys/time.h>


#include "pwlf.h"
#include "filmadastorage.h"
#include "filmadalru.h"
#include "zipf.h"
typedef unsigned long int key_type;  // if (filename== "books" || filename == "wiki_ts", "synthetic", "YCSB")
//typedef double key_type;   //if (filename== "longitudes")

using std::cout;
using std::endl;

#define forceinline inline __attribute__((__always_inline__))
#define ALEX_USE_LZCNT 1

struct AlexCompare {
    template <class T1, class T2>
    bool operator()(const T1& x, const T2& y) const {
        static_assert(
                std::is_arithmetic<T1>::value && std::is_arithmetic<T2>::value,
                "Comparison types must be numeric.");
        return x < y;
    }
};





/*** Helper methods for bitmap ***/

// Extract the rightmost 1 in the binary representation.
// e.g. extract_rightmost_one(010100100) = 000000100
inline uint64_t extract_rightmost_one(uint64_t value) {
    return value & -static_cast<int64_t>(value);
}

// Remove the rightmost 1 in the binary representation.
// e.g. remove_rightmost_one(010100100) = 010100000
inline uint64_t remove_rightmost_one(uint64_t value) {
    return value & (value - 1);
}

// Count the number of 1s in the binary representation.
// e.g. count_ones(010100100) = 3
inline int count_ones(uint64_t value) {
    return static_cast<int>(_mm_popcnt_u64(value));
}

// Get the offset of a bit in a bitmap.
// word_id is the word id of the bit in a bitmap  // by chao (bitmap_pos)
// bit is the word that contains the bit  // by chao (bit_pos)
inline int get_offset(int word_id, uint64_t bit) {
    return (word_id << 6) + count_ones(bit - 1);
}



namespace filminsert{
#define PGM_SUB_EPS(x, epsilon) ((x) <= (epsilon) ? 0 : ((x) - (epsilon)))
//#define PGM_SUB_EPS2(x, epsilon,size) ((x) <= (epsilon) || (x) - (epsilon)>=(size) ? 0 :((x) - (epsilon)) )
#define PGM_ADD_EPS(x, epsilon, size) ((x) + (epsilon) + 2 >= (size) ? (size) : (x) + (epsilon) + 2)

    // record the access_stats, the number of access memory, the number of access disk and the number of cross memory and disk;
    struct access_stats {
        int memnum = 0;
        int disknum =0;
        int crossnum =0;
        int diskpagenum =0;
        int crosspagenum = 0;
        int writenum = 0;
        double querytimeuse =0.0;
        double computetimeuse =0.0;
        double xlookuptime = 0.0;
        double flattenedtime = 0.0;
        double lrutime = 0.0;
        double disktime = 0.0;
        double rdisktime = 0.0;
        double wdisktime = 0.0;
        double sort_piecetime = 0.0;
        unsigned int pagecross = 0;  // 磁盘页 seek 的跨度
        double zipffactor = 0.0;
        string workload = "hotspot";  // "zipf", "random", "hotspot","randomzipf"
        double fill_factor = 0.0;


        void print_stats(){
            struct timeval t1;
            gettimeofday(&t1, NULL);
            cout<<"ID " <<double(t1.tv_sec) << " zipffactor "<< zipffactor  << " memnum " << memnum << " disknum " << disknum << " crossnum "<< crossnum << " diskpagenum ";
            cout<< diskpagenum<<" crosspagenum " << crosspagenum <<" pagecross " << pagecross<<" writenum " << writenum <<" lrutimeuse " <<lrutime <<" disktimeuse " << disktime <<" readdisktimeuse " <<rdisktime
                <<" writedisktimeuse " <<wdisktime << " indexloopuptimeuse " << xlookuptime << " sort_piecetime "  << sort_piecetime << " flattenedtime " << flattenedtime<< " workload " << workload <<"\n";
            ofstream savefile2;
            string performance_file2 = "/home/wamdm/chaohong/clionDir/updatefilm/result/totally_rebuild_reduced_adalru_film_performance.txt";
            savefile2.open(performance_file2,ios::app);
//            savefile << "_ ";
            savefile2<<"ID " <<double(t1.tv_sec) << " zipffactor "<< zipffactor
                     << " memnum " << memnum << " disknum " << disknum << " crossnum "<< crossnum
                    << " diskpagenum ";
            savefile2<< diskpagenum<<" crosspagenum " << crosspagenum <<" pagecross " << pagecross<<" writenum " << writenum
                    <<" lrutimeuse " <<lrutime <<" disktimeuse " <<disktime << " readdisktimeuse " <<rdisktime
                    <<" writedisktimeuse " <<wdisktime << " indexloopuptimeuse " << xlookuptime<< " sort_piecetime "  << sort_piecetime << " flattenedtime " << flattenedtime<< " workload " << workload<<"\n";
            savefile2 << flush;
            savefile2.close();

        }
    };


    // Linear regression model
    template <class T>
    class LinearModel {
    public:
        double a_ = 0;  // slope
        double b_ = 0;  // intercept

        LinearModel() = default;
        LinearModel(double a, double b) : a_(a), b_(b) {}
        explicit LinearModel(const LinearModel& other) : a_(other.a_), b_(other.b_) {}

        void expand(double expansion_factor) {  // not clear to chao
            a_ *= expansion_factor;
            b_ *= expansion_factor;
        }

        inline int predict(T key) const {
//      auto p = static_cast<int>(a_ * static_cast<double>(key) + b_);   // by chao
            return static_cast<int>(a_ * static_cast<double>(key) + b_);
        }

        inline double predict_double(T key) const {
            return a_ * static_cast<double>(key) + b_;
        }
    };


    /**
 * A struct that stores the result of a query to a @ref PGMIndex, that is, a range [@ref lo, @ref hi)
 * centered around an approximate position @ref pos of the sought key.
 */

    inline int cmp(const void *a, const void *b) {
        uint64_t aa = *(uint64_t *) a;
        uint64_t bb = *(uint64_t *) b;
        if (aa < bb)
            return -1;
        if (aa > bb)
            return 1;
        return 0;
    }
    struct ApproxPos {
        int pos; ///< The predicted position of the key.
        int lo;  ///< The lower bound of the range.
        int hi;  ///< The upper bound of the range.
    };

    struct film_stats {
        /// Number of items in the film
        int leaves;
        /// Number of inners in the B+ tree
        vector<int> inners;
        int innernum;
        /// Number of levels
        int numlevel;

        /// Zero initialized
        inline film_stats()
                : leaves(0), numlevel(0) {
        }

    };

    auto key_less_ = AlexCompare();

    // True if a < b
    template <class K>
    forceinline bool key_less(const key_type& a, const K& b) {
        return key_less_(a, b);
    }

// True if a <= b
    template <class K>
    forceinline bool key_lessequal(const key_type& a, const K& b) {
        return !key_less_(b, a);    // b < a , plus ! return false,  b > a , plus ! return true
    }

// True if a > b
    template <class K>
    forceinline bool key_greater(const key_type& a, const K& b) {
        return key_less_(b, a);
    }

// True if a >= b
    template <class K>
    forceinline bool key_greaterequal(const key_type& a, const K& b) {
        return !key_less_(a, b);
    }

// True if a == b
    template <class K>
    forceinline bool key_equal(const key_type& a, const K& b) {
        return !key_less_(a, b) && !key_less_(b, a);
    }




    template <typename key_type,typename value_type, class Compare = AlexCompare,
            class Alloc = std::allocator<key_type>,
            bool allow_duplicates = true>
    class FILMinsert{
    public:

        //    typedef filminsert::FILMinsert< key_type, key_type* > filmadatype;
        typedef filminsert::FILMinsert<key_type, key_type*, AlexCompare, std::allocator<std::pair<key_type, key_type*>>, true> filmadatype;


        Compare key_less_ = Compare();

//        struct Leafpiece;  // 声明leafpiece 结构体

        struct Leafpiece{

            key_type startkey;
            key_type endkey;
            vector<key_type> slotkey;
            vector< bool > locbitmap;   // locbitmap 的作用： location bitmap： check record 是在内存中还是在磁盘中
            vector< void*  > slotdata;   // slotdata 的作用： point to payload or point to evicted table
            vector<bool> delbitmap;   // delbitmap 的作用：check record 是否被删除  初始化全为1

            adalru::localLRU <lruOff_type,key_type* >  intrachain;    //unsigned short int
            Leafpiece* buffer_piece = NULL;    // 用于 存储 任意位置插入的数据
            double slope;
            double intercept;

            Leafpiece() = default;

            Leafpiece(key_type startkey,key_type endkey,double slope, double intercept):
                    startkey(startkey),endkey(endkey),slope(slope),intercept(intercept){};

            explicit Leafpiece(double slope, double intercept):
                    startkey(),endkey(),slope(slope),intercept(intercept){};

            explicit Leafpiece(key_type endkey,double slope, double intercept):
                    startkey(),endkey(endkey),slope(slope),intercept(intercept){};


            inline explicit Leafpiece(const typename internal::insertPWLF<key_type, int>::CanonicalSegment &cs):
                    startkey(cs.get_first_x()),endkey(cs.get_last_x()){
                auto[cs_slope, cs_intercept] = cs.get_floating_point_segment(startkey);
                if (cs_intercept > std::numeric_limits<decltype(intercept)>::max())
                    throw std::overflow_error("Change the type of Segment::intercept to int64");
                slope = cs_slope;
                intercept = cs_intercept;
                slotkey.assign(cs.slotkey.begin(),cs.slotkey.end());
            };

            friend inline bool operator<(const Leafpiece &p, const key_type &key) { return p.key < key; }
            friend inline bool operator<(const key_type &key, const Leafpiece &p) { return key < p.key; }
            friend inline bool operator<(const Leafpiece &p1, const Leafpiece &p2) { return p1.key < p2.key; }

            operator key_type() { return startkey; };

            /**
         * Returns the approximate position of the specified key.
         * @param k the key whose position must be approximated
         * @return the approximate position of the specified key
         */
            forceinline size_t operator()(const key_type &key) const {
                auto pos = int64_t(slope * (key - startkey)) + intercept;
                return pos > 0 ? size_t(pos) : 0ull;
            }

            forceinline void update(typename internal::insertPWLF<key_type, int>::CanonicalSegment &cs){
                startkey = cs.first;
                endkey = cs.last;
                auto[cs_slope, cs_intercept] = cs.get_floating_point_segment(startkey);
                if (cs_intercept > std::numeric_limits<decltype(intercept)>::max())
                    throw std::overflow_error("Change the type of Segment::intercept to int64");
                slope = cs_slope;
                intercept = cs_intercept;
            }
//            m_piece
            std::tuple<bool, unsigned long > delete_key(int pos, bool in_mem_flag){  // 在m_piece 中删除key： // 在内存中还是在磁盘中， 1） 在内存中，更新del_bitmap，更新payload = null,  (key 保持不变， 或者修改为next key， 以适应二分查找)
//                bool xx = delbitmap[pos];
                delbitmap [pos] = true;
                if (in_mem_flag){  // 说明要删除的数据是在内存中
                    auto del_node = (adalru::Node<lruOff_type, key_type *> *) slotdata[pos];  // deleted record, is a pointer, point to a node in local chain
                    auto size = intrachain.remove_node(del_node);
                    slotdata[pos] = NULL;
                    del_node = NULL;
                    return std::make_tuple(size,0);
                }
                else{  // 在磁盘中
//                    auto xx = slotdata[pos];
                    auto short_reduced_evict = (unsigned long) slotdata[pos];

//                    std::cout<<"Lord, You are my refuge from ever to ever" << std::endl;
                    return std::make_tuple(false, short_reduced_evict);  // 返回磁盘的地址
                }

            }

            inline int binary_search_upper_bound(vector<key_type> keys,int l, int r, const key_type key) const {
                while (l < r) {
                    int mid = l + (r - l) / 2;
//                    auto psize = keys.size();
//                    auto piece = keys[mid];// by chao
                    if (key_lessequal((keys[mid]), key)) {
                        l = mid + 1;
                    } else {
                        r = mid;
                    }
                }
                return l;
            }

            // True if a <= b
            template <class K>
            forceinline bool key_lessequal(key_type a, K b) const {
                return !(b<a);    // b < a , plus ! return false,  b > a , plus ! return true
            }

            // True if a > b
            template <class K>
            forceinline bool key_greater(key_type a,  K b)  const {
                return (a>b);
            }

            inline int exponential_search_upper_bound(const vector<key_type> keys,unsigned long m, const key_type key) {
                // Continue doubling the bound until it contains the upper bound. Then use
                // binary search.
                int bound = 1;
                auto keysize = keys.size();
                if (key_lessequal(keys.size(),m))
                    m = keysize -1;
                int l, r;  // will do binary search in range [l, r)
                auto lo =  keys.begin();  //,(--itlevel)->size()
//            auto hi = childpieces.begin();
                auto ckey = (*(lo+m));   // (*leaflo)

                // 读取位置m 的值。internal node gets the start key, leaf nose gets the key stored in m
                if (key_greater((*(lo+m)), key)) {   // if the key in predicted position is greater than the query key
                    int size = m;
                    while (bound < size &&
                           key_greater( (*(lo+m-bound)), key)) {
//                        key_type mkey = (*(lo+m-bound)) ;// chao
//        key_greater(mkey, key);
//                        auto mkeydiff = mkey - key;  //chao
//                    auto mkey2 = key_slots_[m-bound];   //chao
                        bound *= 2;
//                    num_exp_search_iterations_++;
                    }
                    l = m - std::min<int>(bound, size);
                    r = m - bound / 2;
                } else {    // the key in the predicted position is smaller than the query key, by chao
                    int size = keys.size() - m;   // search in the right of the position m, by chao
                    while (bound < size &&
                           key_lessequal((*(lo + m + bound)), key)) {
                        bound *= 2;
//                    num_exp_search_iterations_++;
                    }
                    l = m + bound / 2;
                    r = m + std::min<int>(bound, size);
                }
                return binary_search_upper_bound(keys, l, r, key);
            }


            // the algorithm of merge sort by learned
            template <class lru_type>
            void learnedMergeSort(int Error,lru_type interchain){

                auto curleaf = this;
                interchain->remove(curleaf->startkey);
                while(curleaf->buffer_piece->buffer_piece){
                    curleaf = curleaf->buffer_piece;
//                auto s = this->slotkey.size();  // vect2.insert(vect2.begin(), vect1.begin(), vect1.end());
                    copy( curleaf->slotkey.begin(), curleaf->slotkey.end(), back_inserter(this->slotkey));
                    copy(curleaf->locbitmap.begin(), curleaf->locbitmap.end(), back_inserter(this->locbitmap));
                    copy(curleaf->slotdata.begin(), curleaf->slotdata.end(), back_inserter(this->slotdata));
                    interchain->remove(curleaf->startkey);
                }
                sortPiece * bufferpiece = (sortPiece*) curleaf->buffer_piece;
                vector<key_type> slotkey;
                vector<bool> locbitmap;
                vector<void*> slotdata;
                unsigned long slot = 0;
                // the first key insert use learned sort
                auto key = bufferpiece->slotkey[0];

                //            if (key == 31835864 || key == 31836190)
//                cout << "i need You, my lovely Lord, rebuild 30359458" << endl;
                unsigned long learnedPos = (*this)(key);
                // 根据learned pos 找到大于等于 buffer key 的第一个元素
//            int slotlo = PGM_SUB_EPS2(learnedPos, Error + 1,(this->slotkey.size()-1)) ;  //,(--itlevel)->size()

                int resslot = exponential_search_upper_bound(this->slotkey,learnedPos, key);
                slotkey.insert(slotkey.end(),this->slotkey.begin()+slot,this->slotkey.begin()+resslot);
                slotkey.emplace_back(key);

                slotdata.insert(slotdata.end(),this->slotdata.begin()+slot,this->slotdata.begin()+resslot);
                if (bufferpiece->locbitmap[0])
                    slotdata.emplace_back(this->intrachain.put(0,(key_type*)bufferpiece->slotdata[0]));
                else
                    slotdata.emplace_back(bufferpiece->slotdata[0]);
                locbitmap.insert(locbitmap.end(),this->locbitmap.begin()+slot,this->locbitmap.begin()+resslot);
                locbitmap.emplace_back(bufferpiece->locbitmap[0]);
                slot = resslot;

                for (int i = 1; i< bufferpiece->slope;i ++){
//                auto it = bufferpiece->slotkey.begin(); it !=bufferpiece->slotkey.end(); it++
                    auto key = bufferpiece->slotkey[i];
// 依次比较 buffer 和 this->slotkey 中的元素，哪个小就插入到slotkey 中
                    if (slot < this->slotkey.size()){
                        if (key < this->slotkey[slot]){
                            slotkey.emplace_back(key);
                            if (bufferpiece->locbitmap[i])
                                slotdata.emplace_back(this->intrachain.put(0,(key_type*)bufferpiece->slotdata[i]));
                            else
                                slotdata.emplace_back(bufferpiece->slotdata[i]);
                            locbitmap.emplace_back(bufferpiece->locbitmap[i]);
                        } else{
                            i --;
                            slotkey.emplace_back(this->slotkey[slot]);
                            slotdata.emplace_back(this->slotdata[slot]);
                            locbitmap.emplace_back(this->locbitmap[slot++]);
                        }
                    }
                    else{
                        slotkey.insert(slotkey.end(),bufferpiece->slotkey.begin()+i,bufferpiece->slotkey.begin()+bufferpiece->slotkey.size());
                        locbitmap.insert(locbitmap.end(),bufferpiece->locbitmap.begin()+i,bufferpiece->locbitmap.begin()+bufferpiece->locbitmap.size());
                        for (;i<bufferpiece->slope;i ++){
                            if (bufferpiece->locbitmap[i])
                                slotdata.emplace_back(this->intrachain.put(0,(key_type*)bufferpiece->slotdata[i]));
                            else
                                slotdata.emplace_back(bufferpiece->slotdata[i]);
                        }
                        break;
                    }
                    // i 为元素对应的下标
                }
//            auto nowkey = this->slotkey[slot];
                slotkey.insert(slotkey.end(),this->slotkey.begin()+slot,this->slotkey.end());
                slotdata.insert(slotdata.end(),this->slotdata.begin()+slot,this->slotdata.end());
                locbitmap.insert(locbitmap.end(),this->locbitmap.begin()+slot,this->locbitmap.end());

                this->slotkey = slotkey;
                this->slotdata = slotdata;
                this->locbitmap = locbitmap;

                bufferpiece->slope = 0;
                bufferpiece->intercept = 0;
//            cout << "i need You, my Lord" << endl;
//            要将bufferpiece 从interchain 也就是global LRU中删除
                interchain->remove(bufferpiece->slotkey[0]);  // 要将bufferpiece 从interchain 也就是global LRU中删除
            }

            template <class lru_type>   // 这个需要考虑，如果已经被删除了的key, 就不参与 retrain 了
            void mlearnedMergeSort(int Error,lru_type interchain,std::vector<key_type>& retrainkeys, std::vector<void*>& retraindata,std::vector<bool>& retrainflag){

                auto curleaf = this;
                interchain->remove(curleaf->startkey);
                while(curleaf->buffer_piece->buffer_piece){
                    curleaf = curleaf->buffer_piece;
//                auto s = this->slotkey.size();  // vect2.insert(vect2.begin(), vect1.begin(), vect1.end());
                    copy( curleaf->slotkey.begin(), curleaf->slotkey.end(), back_inserter(this->slotkey));
                    copy(curleaf->locbitmap.begin(), curleaf->locbitmap.end(), back_inserter(this->locbitmap));
                    copy(curleaf->delbitmap.begin(), curleaf->delbitmap.end(), back_inserter(this->delbitmap));
                    copy(curleaf->slotdata.begin(), curleaf->slotdata.end(), back_inserter(this->slotdata));
                    interchain->remove(curleaf->startkey);
                }
                sortPiece * bufferpiece = (sortPiece*) curleaf->buffer_piece;
//                vector<key_type> retrainkeys;
//                vector<bool> retrainflag;
//                vector<void*> retraindata;
                unsigned long slot = 0;
                // the first key insert use learned sort
                auto key = bufferpiece->ALEX_DATA_NODE_KEY_AT(0);

//                if (key == 31835864 || key == 31836190)
//                    cout << "i need You, my lovely Lord, rebuild 30359458" << endl;
                unsigned long learnedPos = (*this)(key);
                // 根据learned pos 找到大于等于 buffer key 的第一个元素
//            int slotlo = PGM_SUB_EPS2(learnedPos, Error + 1,(this->slotkey.size()-1)) ;  //,(--itlevel)->size()

                int resslot = exponential_search_upper_bound(this->slotkey,learnedPos, key);
                retrainkeys.insert(retrainkeys.end(),this->slotkey.begin()+slot,this->slotkey.begin()+resslot);
                retrainkeys.emplace_back(key);

                retraindata.insert(retraindata.end(),this->slotdata.begin()+slot,this->slotdata.begin()+resslot);
                if (bufferpiece->ALEX_DATA_NODE_FLAG_AT(0))
                    retraindata.emplace_back(this->intrachain.put(0,(key_type*)bufferpiece->ALEX_DATA_NODE_PAYLOAD_AT(0)));
                else
                    retraindata.emplace_back(bufferpiece->ALEX_DATA_NODE_PAYLOAD_AT(0));

                retrainflag.insert(retrainflag.end(),this->locbitmap.begin()+slot,this->locbitmap.begin()+resslot);
                retrainflag.emplace_back(bufferpiece->ALEX_DATA_NODE_FLAG_AT(0));
                slot = resslot;

                for (int i = 1; i< bufferpiece->num_keys_;i ++){
                    auto key = bufferpiece->key_slots_[i];
// 依次比较 buffer 和 this->slotkey 中的元素，哪个小就插入到slotkey 中
                    if (slot < this->slotkey.size()) {
                        if (key < this->slotkey[slot]) {
                            retrainkeys.emplace_back(key);
                            if (bufferpiece->pointer_slots_[i].first)
                                retraindata.emplace_back(
                                        this->intrachain.put(0, (key_type *) bufferpiece->pointer_slots_[i].second));
                            else
                                retraindata.emplace_back(bufferpiece->pointer_slots_[i].second);
                            retrainflag.emplace_back(bufferpiece->pointer_slots_[i].first);
                        } else {
                            i--;
                            retrainkeys.emplace_back(this->slotkey[slot]);
                            retraindata.emplace_back(this->slotdata[slot]);
                            retrainflag.emplace_back(this->locbitmap[slot++]);
                        }
                    } else {
                        retrainkeys.insert(retrainkeys.end(), bufferpiece->key_slots_ + i,
                                       bufferpiece->key_slots_ + bufferpiece->num_keys_);
                        for (; i < bufferpiece->num_keys_; i++) {
                            retrainflag.emplace_back(bufferpiece->pointer_slots_[i].first);
                            if (bufferpiece->pointer_slots_[i].first)
                                retraindata.emplace_back(
                                        this->intrachain.put(0, (key_type *) bufferpiece->pointer_slots_[i].second));
                            else
                                retraindata.emplace_back(bufferpiece->pointer_slots_[i].second);
                        }
                        break;
                    }
                    // i 为元素对应的下标
                }
//            auto nowkey = this->slotkey[slot];
                retrainkeys.insert(retrainkeys.end(),this->slotkey.begin()+slot,this->slotkey.end());
                retraindata.insert(retraindata.end(),this->slotdata.begin()+slot,this->slotdata.end());
                retrainflag.insert(retrainflag.end(),this->locbitmap.begin()+slot,this->locbitmap.end());

//                this->slotkey = slotkey;
//                this->slotdata = slotdata;
//                this->locbitmap = locbitmap;

//            cout << "i need You, my Lord" << endl;
//            要将bufferpiece 从interchain 也就是global LRU中删除
                interchain->remove(bufferpiece->ALEX_DATA_NODE_KEY_AT(0));  // 要将bufferpiece 从interchain 也就是global LRU中删除
                vector<key_type>().swap(this->slotkey);
                vector<bool>().swap(this->locbitmap);
                vector<void*>().swap(this->slotdata);
                /*
// 需要保留的代码，如果是根据阈值来retrain，就需要check_exists. 满了以后再retrain，就不需要判断
for (int i = 0; i< bufferpiece->data_capacity_;i ++){
if (bufferpiece->check_exists(i)){
auto key = bufferpiece->key_slots_[i];
// 依次比较 buffer 和 this->slotkey 中的元素，哪个小就插入到slotkey 中
if (slot < this->slotkey.size()){
if (key < this->slotkey[slot]){
slotkey.emplace_back(key);
if (bufferpiece->pointer_slots_[i].first)
slotdata.emplace_back(this->intrachain.put(0,(key_type*)bufferpiece->pointer_slots_[i].second));
else
slotdata.emplace_back(bufferpiece->pointer_slots_[i].second);
locbitmap.emplace_back(bufferpiece->pointer_slots_[i].first);
} else{
i --;
slotkey.emplace_back(this->slotkey[slot]);
slotdata.emplace_back(this->slotdata[slot]);
locbitmap.emplace_back(this->locbitmap[slot++]);
}
}
else{
slotkey.insert(slotkey.end(),bufferpiece->key_slots_+i,bufferpiece->key_slots_+bufferpiece->num_keys_);
//                        locbitmap.insert(locbitmap.end(),bufferpiece->locbitmap.begin()+i,bufferpiece->locbitmap.begin()+bufferpiece->locbitmap.size());
for (;i<bufferpiece->num_keys_;i ++){
if (bufferpiece->check_exists(i)){
locbitmap.emplace_back(bufferpiece->pointer_slots_[i].first);
if (bufferpiece->pointer_slots_[i].first)
slotdata.emplace_back(this->intrachain.put(0,(key_type*)bufferpiece->pointer_slots_[i].second));
else
slotdata.emplace_back(bufferpiece->pointer_slots_[i].second);
}
}
break;
}
}
// i 为元素对应的下标
}
 */

            }



            // 判断key 是否在当前leaf 中
            inline pair<bool,unsigned long> is_in_leaf(const key_type &key,const int error){
                pair<bool,unsigned long> is_in(false,0);
                if (key<=this->endkey ){
                    auto pos = (*this)(key);

                    int slotlo = PGM_SUB_EPS(pos, error + 1) ;  //,(--itlevel)->size()
                    int slothi = PGM_ADD_EPS(pos, error,(this->slotkey.size()-1));
                    for (; slotlo<=slothi;slotlo++){
//                    if (slotlo < 0)
//                    {
//                        cout<< "my Lord, please help, i trust in You!" << endl;
//                    }
                        if (key < this->slotkey[slotlo] || key >=this->startkey)
                        {
                            is_in.first = true;
                            is_in.second = slotlo;
                            break;
                        }
                        else
                        {
                            if (key == this->endkey && key<this->startkey){
                                is_in.first = true;
                                is_in.second = slotlo;
                                break;
                            }
//                        auto aaa = this->slotkey[slotlo];
//                        cout<< "my Lord, please help, i trust in You!" << endl;
                        }
//                    else{
//                        // check in sort piece
//                        cout << "i need You, my Lord" << endl;
////                        auto sortpos = (*sort_list)(key);
////                        if (sort_list->slotkey[sortpos] == key){
////                            is_in.first = true;
////                            is_in.second = slotlo;
////                            break;
////                        }
//                    }
                    }

                    return is_in;
                }
                else
                    return is_in;
            }

            // 找到key 在当前leaf 中的upper position, that is the last key that larger or equal to the key
            inline pair<bool,unsigned long> find_upper_in_leaf(const key_type &key,const int error){
                pair<bool,unsigned long> is_in(false,this->slotkey.size()-1);
                if (key<=this->endkey ){
                    if (key < startkey ){
                        is_in.first = true;
                        return is_in;
//                    cout << "Jesus, i trust in You" << endl;
                    }
                    auto pos = (*this)(key);

                    unsigned int slotlo = PGM_SUB_EPS(pos, error + 1) ;  //,(--itlevel)->size()
                    unsigned int slothi = PGM_ADD_EPS(pos, error,(this->slotkey.size()-1));
                    for (; slotlo<=slothi;slotlo++){

                        if (key <= this->slotkey[slotlo])
                        {
                            is_in.first = true;
//                            auto upperkey = this->slotkey[slotlo];
                            is_in.second = slotlo;
                            break;
                        }

                    }

                    return is_in;
                }
                else
                    return is_in;
            }

        };


        struct sortPiece:Leafpiece{

            typedef std::pair<key_type, value_type> V;
//            typedef sortPiece<key_type, value_type, Compare, Alloc, allow_duplicates> self_type;
            typedef typename Alloc::template rebind<sortPiece>::other alloc_type;
            typedef typename Alloc::template rebind<key_type>::other key_alloc_type;
            typedef typename Alloc::template rebind<value_type>::other payload_alloc_type;
            typedef typename Alloc::template rebind<std::pair<bool,void*>>::other pointer_alloc_type;  //by chao
            typedef typename Alloc::template rebind<V>::other value_alloc_type;
            typedef typename Alloc::template rebind<uint64_t>::other bitmap_alloc_type;
            static constexpr key_type kEndSentinel_ = std::numeric_limits<key_type>::max();
            const Compare& key_less_ = Compare();
            const Alloc& allocator_ =  Alloc();
//            Compare key_less_ = Compare();
//            Alloc allocator_ = Alloc();

//            // Forward declaration
//            template <typename node_type = self_type, typename payload_return_type = P,
//                    typename value_return_type = V>
//            class Iterator;
//            typedef Iterator<> iterator_type;
//            typedef Iterator<const self_type, const P, const V> const_iterator_type;

//            LinearModel<key_type>* model_;
            uint64_t* bitmap_ = nullptr;
            int bitmap_size_ = 0;  // number of int64_t in bitmap

            int data_capacity_ = 0;  // size of key/data_slots array
            int num_keys_ = 0;  // number of filled key/data slots (as opposed to gaps)
            int num_shifts_ = 0;
            vector<key_type> slotkey;
//        vector< pair<bool, adalru::Node<key_type, std::vector<key_type> >*  > >slotdata;
            vector< bool> locbitmap;   //adalru::Node<key_type, std::vector<key_type> >*sss
            vector< void* > slotdata;   //adalru::Node<key_type, std::vector<key_type> >*sss
//            adalru::localLRU <lruOff_type,key_type* >  intrachain;    //unsigned short int
            key_type* key_slots_ = nullptr;  // holds keys
//        P* payload_slots_ =
//                nullptr;  // holds payloads, must be same size as key_slots//
            std::pair<bool, void*>* pointer_slots_ = nullptr;
            ~sortPiece() {
                delete key_slots_;// delete
                delete pointer_slots_;
                delete bitmap_;
            }
            sortPiece(int capacity):data_capacity_(capacity){
//                key_less_(other.key_less_),
//                        allocator_(other.allocator_),
                bitmap_size_ = static_cast<size_t>(std::ceil(data_capacity_ / 64.));
                bitmap_ = new (bitmap_allocator().allocate(bitmap_size_))
                        uint64_t[bitmap_size_]();  // initialize to all false
                key_slots_ =
                        new (key_allocator().allocate(data_capacity_)) key_type[data_capacity_];
//                memset(key_slots_, init_key,sizeof(key_slots_) );
                for (int i = 0; i < data_capacity_; i++) {  // by chao, fill the gaps with the closest key to the right of the gap, which helps maintain exponential search performance
                    ALEX_DATA_NODE_KEY_AT(i) = kEndSentinel_;
                }
//            payload_slots_ =
//                    new (payload_allocator().allocate(data_capacity_)) P[data_capacity_];
                pointer_slots_ =
                        new (pointer_allocator().allocate(data_capacity_)) std::pair<bool, void*>[data_capacity_];
            }

            key_alloc_type key_allocator() { return key_alloc_type(allocator_); }

            payload_alloc_type payload_allocator() {
                return payload_alloc_type(allocator_);
            }
            pointer_alloc_type pointer_allocator() {
                return pointer_alloc_type(allocator_);
            }

            value_alloc_type value_allocator() { return value_alloc_type(allocator_); }

            bitmap_alloc_type bitmap_allocator() { return bitmap_alloc_type(allocator_); }



            friend inline bool operator<(const sortPiece &p, const key_type &key) { return p.key < key; }
            friend inline bool operator<(const key_type &key, const sortPiece &p) { return key < p.key; }
            friend inline bool operator<(const sortPiece &p1, const sortPiece &p2) { return p1.key < p2.key; }


            /**
         * Returns the approximate position of the specified key.
         * @param k the key whose position must be approximated
         * @return the approximate position of the specified key
         */
            forceinline int operator()(const key_type key)  const{
                int begin = 0, end = slotkey.size()-1,mid;  // find last Equal or small, 找到最后一个大于或等于 key 的
                while(begin <= end){
                    mid = (end + begin) / 2;
                    if(slotkey[mid] >= key) {
                        end = mid-1;
                    } else
                        begin = mid +1;
                }
                return end+1;
            }


            forceinline pair<bool,bool>  insert(key_type key, key_type* payload) {
                // the first bool flag is used to indicate rebuild,the second is used to start key changing
                int begin = 0, end = slotkey.size()-1,mid;  // find last Equal or small, 找到最后一个大于或等于 key 的
                while(begin <= end){
                    mid = (end + begin) / 2;
                    if(this->slotkey[mid] >= key) {
                        end = mid-1;
                    } else
                        begin = mid +1;
                }

//                if (end < 0 && this->slotkey.size() >0)// 所有 在 end 之后的，都需要加1
                this->slotkey.insert(slotkey.begin()+end+1,key);
                this->slotdata.insert( slotdata.begin()+end+1,payload);
                this->locbitmap.insert(locbitmap.begin()+end+1,true);
                if (++this->intercept == this->slope ){
                    // rebuild this leaf piece, merge the slotkey and buffer piece
                    if (end < 0 && this->slotkey.size() >1){
//                        cout << "i need You, my Lord, please come, rebuild this piece, and change the startkey" << endl;
                        return pair<bool,bool> (true,true);
                    }
                    else{
//                        cout << "i need You, my Lord, please come, rebuild this piece" << endl;
                        return pair<bool,bool> (true,false);
                    }
                }
                else{
                    if (this->intercept != this->num_keys_)
                        COUT_THIS("Jesus, please come");
                    if (end < 0 && this->slotkey.size() > 1){
//                        cout<< "i need You, my Lord" << endl;
                        return pair<bool,bool> (false,true);
                    }// 这表示start key 发生了变化
                    else{
                        return pair<bool,bool> (false,false);
                    }

                }

            }


            // Predicts the position of a key using the model
            inline int predict_position(const key_type& key) const {
                int position = this->model_.predict(key);
                position = std::max<int>(std::min<int>(position, data_capacity_ - 1), 0);
                return position;
            }

            // Check whether the position corresponds to a key (as opposed to a gap)
            inline bool check_exists(int pos) const {
//                assert(pos >= 0 && pos < data_capacity_);
                int bitmap_pos = pos >> 6;  // divide by 64, by chao
                int bit_pos = pos - (bitmap_pos << 6);
                return static_cast<bool>(bitmap_[bitmap_pos] & (1ULL << bit_pos));
            }

            // Mark the entry for position in the bitmap
            inline void set_bit(int pos) {
                assert(pos >= 0 && pos < data_capacity_);
                int bitmap_pos = pos >> 6;    // by chao, bitmap_pos represents the position of uint64
                int bit_pos = pos - (bitmap_pos << 6);     // by chao, bit_pos represents the position in uint64

                bitmap_[bitmap_pos] |= (1ULL << bit_pos); // by chao, 将unsigned long long 1 (64-bit 1), left shift bit_pos, 这表示，只将 bit_pos 位置为1， 其他位，采用 | 或操作，保持不变
            }

            // Mark the entry for position in the bitmap
            inline void set_bit(uint64_t bitmap[], int pos) {
                int bitmap_pos = pos >> 6;
                int bit_pos = pos - (bitmap_pos << 6);
                bitmap[bitmap_pos] |= (1ULL << bit_pos);
            }

            // Unmark the entry for position in the bitmap
            inline void unset_bit(int pos) {
                assert(pos >= 0 && pos < data_capacity_);
                int bitmap_pos = pos >> 6;
                int bit_pos = pos - (bitmap_pos << 6);
                bitmap_[bitmap_pos] &= ~(1ULL << bit_pos);
            }

            // find the first valid record position in the bitmap
            inline int find_first_1() {
                for (int i = 0; i < bitmap_size_; i++){
                    uint64_t slice = bitmap_[i];

//                    auto bit_pos = binarysearch1(slice);
//                    auto bit0 = bsr_asm(0);
//                    auto bit1 = bsr_asm(1);
//                    auto bit2 = bsr_asm(2);
                    if (slice){
//                        auto bit_pos1 = bsr_asm(slice);
                        auto bit_pos = bsf_asm(slice);
                        int first_1_pos  = (i)*64 + bit_pos;
//                        if (bit_pos1 != bit_pos ){
//                            COUT_THIS("thank You, my lovely Lord");
//                        }
                        return first_1_pos; //1000，0000，0000，0000，0000，0000，0000
                    }
                }
                return data_capacity_;


            }

            // find the first valid record position in the bitmap
            inline int find_first_1_incorrect(int pos) {  // 如果当前slice 存在多个valid，那么就会出错，尤其当valid 小于 pos 时
                int i = pos >> 6;
                for (; i < bitmap_size_; i++){
                    uint64_t slice = bitmap_[i];

//                    auto bit_pos = binarysearch1(slice);
//                    auto bit0 = bsr_asm(0);
//                    auto bit1 = bsr_asm(1);
//                    auto bit2 = bsr_asm(2);
                    if (slice){
//                        auto bit_pos1 = bsr_asm(slice);
                        auto bit_pos = bsf_asm(slice);
                        int first_1_pos  = (i)*64 + bit_pos;
//                        if (bit_pos1 != bit_pos ){
//                            COUT_THIS("thank You, my lovely Lord");
//                        }
                        if (first_1_pos >= pos)
                            return first_1_pos; //1000，0000，0000，0000，0000，0000，0000
                    }
                }
                return data_capacity_;


            }

            // find the next valid record position in the bitmap
            inline int find_first_1(int pos) {
                int i = pos >> 6;
                uint64_t slice = bitmap_[i];
                auto bit_pos = bsr_asm(slice);
                int bit_pos0 = pos - (i << 6);
                slice = slice&(slice-1);   // 移除最右侧的1，然后再找第一个1
                // 先判断最右侧的1 是否是pos，如果是，移除后，再找第一个1，如果不是，移除，直到把pos位的1 移除掉
                if (slice){
                    auto bit_pos = bsr_asm(slice);
                    while (slice && bit_pos <= bit_pos0 ){
                        bit_pos = bsr_asm(slice);
                        slice = slice&(slice-1);
                    }
                    if (bit_pos > bit_pos0) {
                        int first_1_pos = (i) * 64 + bit_pos;
                        return first_1_pos;
                    }
                    else{
                        i ++;
                        // 对slice 中的 大于 bit_pos 的所有bits 逐个判断是否为1
                        for (; i < bitmap_size_; i++){
                            uint64_t slice = bitmap_[i];

//                    auto bit_pos = binarysearch1(slice);
//                    auto bit0 = bsr_asm(0);
//                    auto bit1 = bsr_asm(1);
//                    auto bit2 = bsr_asm(2);
                            if (slice){
//                        auto bit_pos1 = bsr_asm(slice);
                                auto bit_pos = bsf_asm(slice);
                                int first_1_pos  = (i)*64 + bit_pos;
//                        if (bit_pos1 != bit_pos ){
//                            COUT_THIS("thank You, my lovely Lord");
//                        }

                                return first_1_pos; //1000，0000，0000，0000，0000，0000，0000
                            }
                        }
                    }
                }
                else{
                    i ++;
                    // 对slice 中的 大于 bit_pos 的所有bits 逐个判断是否为1
                    for (; i < bitmap_size_; i++){
                        uint64_t slice = bitmap_[i];

//                    auto bit_pos = binarysearch1(slice);
//                    auto bit0 = bsr_asm(0);
//                    auto bit1 = bsr_asm(1);
//                    auto bit2 = bsr_asm(2);
                        if (slice){
//                        auto bit_pos1 = bsr_asm(slice);
                            auto bit_pos = bsf_asm(slice);
                            int first_1_pos  = (i)*64 + bit_pos;
//                        if (bit_pos1 != bit_pos ){
//                            COUT_THIS("thank You, my lovely Lord");
//                        }

                            return first_1_pos; //1000，0000，0000，0000，0000，0000，0000
                        }
                    }
                }


                return data_capacity_;


            }

            // find the last valid record position in the bitmap
            inline int find_last_1() {
                for (int i = bitmap_size_-1; i >= 0; i--){
                    uint64_t slice = bitmap_[i];

//                    auto bit_pos = binarysearch1(slice);
//                    auto bit0 = bsr_asm(0);
//                    auto bit1 = bsr_asm(1);
//                    auto bit2 = bsr_asm(2);
                    if (slice){
                        auto bit_pos = bsr_asm(slice);
                        auto bit_pos1 = bsf_asm(slice);
                        int last_1_pos  = i*64 + bit_pos;
//                        if (bit_pos1 != bit_pos ){
//                            COUT_THIS("thank You, my lovely Lord");
//                        }
                        return last_1_pos; //1000，0000，0000，0000，0000，0000，0000
                    }
                }
                return 0;
            }



            inline int bsr_asm (uint64_t w)  // from the low to high (right to left)
            {
                   int x1, x2;
                   asm ("bsr %1,%0\n" "jnz 1f\n" "bsr %0,%0\n" "subl $32,%0\n"
                      "1: addl $32,%0\n": "=&q" (x1), "=&q" (x2):"1" ((int) (w >> 32)),
                              "0" ((int) w));
                  return x1;
            }

             inline int bsf_asm (uint64_t w)
             {
                   int x1, x2;
                   asm ("bsf %0,%0\n" "jnz 1f\n" "bsf %1,%0\n" "jz 1f\n" "addl $32,%0\n"
                      "1:": "=&q" (x1), "=&q" (x2):"1" ((int) (w >> 32)),
                      "0" ((int) w));
                   return x1;
             }

            // Searches for the last non-gap position equal to key
            // If no positions equal to key, returns -1
            int find_key(const key_type& key,int predicted_pos) {
                predicted_pos = std::min<int>(predicted_pos, data_capacity_ - 1);
//                int predicted_pos = predict_position(key);

                // The last key slot with a certain value is guaranteed to be a real key
                // (instead of a gap)
//                if (key == ALEX_DATA_NODE_KEY_AT(predicted_pos)){  //
//                    if (!key_equal(ALEX_DATA_NODE_KEY_AT(predicted_pos), key)) {
//                        return -1;
//                    } else {
//                        return predicted_pos;
//                    }
//                }
                int pos = exponential_search_upper_bound(predicted_pos, key) - 1;

                if (pos < 0 || !key_equal(ALEX_DATA_NODE_KEY_AT(pos), key)) {
                    return -1;
                } else {
                    return pos;
                }
            }

            // b_piece: 在指定的位置删除key，及对应的record
            unsigned long delete_key(int pos,bool in_mem_flag){
//                predicted_pos = std::min<int>(predicted_pos, data_capacity_ - 1);
//                int pos = exponential_search_upper_bound(predicted_pos, key) - 1;
                // 找到pos 之前，最近的一个有效位置
                auto del_key = this->ALEX_DATA_NODE_KEY_AT(pos);
                auto del_pos = pos-1;
                auto next_pos = pos+1;
                auto next_key = this->ALEX_DATA_NODE_KEY_AT(next_pos);

//                auto xxx = ALEX_DATA_NODE_PAYLOAD_AT(pos);
                while (del_pos >= 0 && ! this->check_exists(del_pos)) {
                    del_pos--;
//                    xxx = ALEX_DATA_NODE_PAYLOAD_AT(del_pos);
//                    std::cout<<"Lord, You are my refuge from ever to ever" << std::endl;
                }
                // std::fill 是左闭右开
                std::fill(&this->ALEX_DATA_NODE_KEY_AT(del_pos + 1),
                          &this->ALEX_DATA_NODE_KEY_AT(next_pos),
                          next_key);

                // 如果在内存中，首先，释放payload所在的指针，并置为null
                if (in_mem_flag){
                    ALEX_DATA_NODE_PAYLOAD_AT(pos) = NULL;
                    this->ALEX_DATA_NODE_FLAG_AT(pos) = false;
//                    std::cout<<"Lord, You are my refuge from ever to ever" << std::endl;
                    return 0;
                }
                else{  // 要删除的record 是在磁盘中， 首先需要读取磁盘中的地址；修改page information
//                    auto xx = this->ALEX_DATA_NODE_PAYLOAD_AT(pos);
                    auto short_reduced_evict = (unsigned long) this->ALEX_DATA_NODE_PAYLOAD_AT(pos);
                    ALEX_DATA_NODE_PAYLOAD_AT(pos) = NULL;
//                    std::cout<<"Lord, You are my refuge from ever to ever" << std::endl;
                    return short_reduced_evict;  // 返回磁盘的地址
                }

//                while (!del_flag ){
//                    del_key = this->ALEX_DATA_NODE_KEY_AT(pos);
//                    this->ALEX_DATA_NODE_KEY_AT(pos) = next_key;
//                    pos--;
//                    del_flag = this->ALEX_DATA_NODE_FLAG_AT(pos);
//                    std::cout<< "Lord, please come, chao needs You" << std::endl;
//                }
            }

            int find_key_lower(const key_type& key,int predicted_pos) {
                predicted_pos = std::max<int>(std::min<int>(predicted_pos, data_capacity_ - 1), 0);
//                int predicted_pos = predict_position(key);

                // The last key slot with a certain value is guaranteed to be a real key
                // (instead of a gap)
//                if (key == ALEX_DATA_NODE_KEY_AT(predicted_pos)){  //
//                    if (!key_equal(ALEX_DATA_NODE_KEY_AT(predicted_pos), key)) {
//                        return -1;
//                    } else {
//                        return predicted_pos;
//                    }
//                }
                int pos = exponential_search_upper_bound(predicted_pos, key) - 1;
                // check 一下 pos 是否valid，如果不valid，找到valid 的下一个位置
                return pos;

            }

            int find_key_upper(const key_type& key,int predicted_pos) {
                predicted_pos = std::max<int>(std::min<int>(predicted_pos, data_capacity_ - 1), 0);
//                int predicted_pos = predict_position(key);

                // The last key slot with a certain value is guaranteed to be a real key
                // (instead of a gap)
//                if (key == ALEX_DATA_NODE_KEY_AT(predicted_pos)){  //
//                    if (!key_equal(ALEX_DATA_NODE_KEY_AT(predicted_pos), key)) {
//                        return -1;
//                    } else {
//                        return predicted_pos;
//                    }
//                }
                int pos = exponential_search_upper_bound(predicted_pos, key);
//                auto mkey = ALEX_DATA_NODE_KEY_AT(pos);
//                while (ALEX_DATA_NODE_KEY_AT(pos)< key){
//                    mkey = ALEX_DATA_NODE_KEY_AT(pos);
//                    pos++;
//                }

//                mkey = ALEX_DATA_NODE_KEY_AT(pos);
                return pos;

            }

            int find_lastkey(const key_type& key,int predicted_pos) {
                predicted_pos = std::max<int>(std::min<int>(predicted_pos, data_capacity_ - 1), 0);
//                int predicted_pos = predict_position(key);

                // The last key slot with a certain value is guaranteed to be a real key
                // (instead of a gap)
                while(!key_equal(ALEX_DATA_NODE_KEY_AT(predicted_pos), key) || !check_exists(predicted_pos)){
                    predicted_pos++;
                }
                if (!key_equal(ALEX_DATA_NODE_KEY_AT(predicted_pos), key))
                    return -1;
                return predicted_pos;
            }

            // Searches for the first non-gap position no less than key
            // Returns position in range [0, data_capacity]
            // Compare with lower_bound()
            int find_lower(const key_type& key) {

                int predicted_pos = predict_position(key);

                int pos = exponential_search_lower_bound(predicted_pos, key);
                return get_next_filled_position(pos, false);
            }

            // Searches for the first non-gap position greater than key
            // Returns position in range [0, data_capacity]
            // Compare with upper_bound()
            int find_upper(const key_type& key) {

                int predicted_pos = predict_position(key);

                int pos = exponential_search_upper_bound(predicted_pos, key);
                return get_next_filled_position(pos, false);
            }

            // Finds position to insert a key.
            // First returned value takes prediction into account.
            // Second returned value is first valid position (i.e., upper_bound of key).
            // If there are duplicate keys, the insert position will be to the right of
            // all existing keys of the same value.
            std::pair<int, int> find_insert_position(const key_type& key, int predicted_pos) {
//                int predicted_pos = predict_position(key);  // first use model to get prediction

//                auto findk = ALEX_DATA_NODE_KEY_AT(predicted_pos);
//                auto findp = ALEX_DATA_NODE_PAYLOAD_AT(predicted_pos);
//                int insert_pos = 0;
//                if (ALEX_DATA_NODE_KEY_AT(predicted_pos) > key){
//                    COUT_THIS("please come, my Lord!");
//                }
                // insert to the right of duplicate keys
                int pos = exponential_search_upper_bound(predicted_pos, key);


                if (predicted_pos <= pos || check_exists(pos)) {
                    return {pos, pos};
                } else {
                    // Place inserted key as close as possible to the predicted position while
                    // maintaining correctness
                    return {std::min(predicted_pos, get_next_filled_position(pos, true) - 1),
                            pos};
                }
            }

            // Starting from a position, return the first position that is not a gap
            // If no more filled positions, will return data_capacity
            // If exclusive is true, output is at least (pos + 1)
            // If exclusive is false, output can be pos itself
            int get_next_filled_position(int pos, bool exclusive) const {
                if (exclusive) {
                    pos++;
                    if (pos == data_capacity_) {
                        return data_capacity_;
                    }
                }

                int curBitmapIdx = pos >> 6;
                uint64_t curBitmapData = bitmap_[curBitmapIdx];

                // Zero out extra bits
                int bit_pos = pos - (curBitmapIdx << 6);
                curBitmapData &= ~((1ULL << (bit_pos)) - 1);

                while (curBitmapData == 0) {
                    curBitmapIdx++;
                    if (curBitmapIdx >= bitmap_size_) {
                        return data_capacity_;
                    }
                    curBitmapData = bitmap_[curBitmapIdx];
                }
                uint64_t bit = extract_rightmost_one(curBitmapData);
                return get_offset(curBitmapIdx, bit);
            }

            // Searches for the first position greater than key
            // This could be the position for a gap (i.e., its bit in the bitmap is 0)
            // Returns position in range [0, data_capacity]
            // Compare with find_upper()
            template <class K>
            int upper_bound(const K& key) {

                int position = predict_position(key);
                return exponential_search_upper_bound(position, key);
            }

            // Searches for the first position greater than key, starting from position m
            // Returns position in range [0, data_capacity]
            template <class K>
            inline int exponential_search_upper_bound(int m, const K& key) {
                // Continue doubling the bound until it contains the upper bound. Then use
                // binary search.
                int bound = 1;
                int l, r;  // will do binary search in range [l, r)
                if (key_greater(ALEX_DATA_NODE_KEY_AT(m), key)) {   // if the key in predicted position is greater than the query key
                    int size = m;
                    while (bound < size &&
                           key_greater(ALEX_DATA_NODE_KEY_AT(m - bound), key)) {
//                        key_type mkey = ALEX_DATA_NODE_KEY_AT(m - bound)  ;// chao
//        key_greater(mkey, key);  // 12061975546076
//                        auto mkeydiff = mkey - key;  //chao
//                        auto mkey2 = key_slots_[m-bound];   //chao
                        bound *= 2;
                    }
                    l = m - std::min<int>(bound, size);
                    r = m - bound / 2;
                } else {    // the key in the predicted position is smaller than the query key, by chao
                    int size = data_capacity_ - m;   // search in the right of the position m, by chao
                    while (bound < size &&
                           key_lessequal(ALEX_DATA_NODE_KEY_AT(m + bound), key)) {
                        bound *= 2;
                    }
                    l = m + bound / 2;
                    r = m + std::min<int>(bound, size);
                }
                return binary_search_upper_bound(l, r, key);
            }


            // Searches for the first position greater than key in range [l, r)
            // https://stackoverflow.com/questions/6443569/implementation-of-c-lower-bound
            // Returns position in range [l, r]
            template <class K>
            inline int binary_search_upper_bound(int l, int r, const K& key) const {
                while (l < r) {
                    int mid = l + (r - l) / 2;
                    if (key_lessequal(ALEX_DATA_NODE_KEY_AT(mid), key)) {
                        l = mid + 1;
                    } else {
                        r = mid;
                    }
                }
                return l;
            }

            // Searches for the first position no less than key
            // This could be the position for a gap (i.e., its bit in the bitmap is 0)
            // Returns position in range [0, data_capacity]
            // Compare with find_lower()
            template <class K>
            int lower_bound(const K& key) {

                int position = predict_position(key);
                return exponential_search_lower_bound(position, key);
            }

            // Searches for the first position no less than key, starting from position m
            // Returns position in range [0, data_capacity]
            template <class K>
            inline int exponential_search_lower_bound(int m, const K& key) {
                // Continue doubling the bound until it contains the lower bound. Then use
                // binary search.
                int bound = 1;
                int l, r;  // will do binary search in range [l, r)
                if (key_greaterequal(ALEX_DATA_NODE_KEY_AT(m), key)) {
                    int size = m;
                    while (bound < size &&
                           key_greaterequal(ALEX_DATA_NODE_KEY_AT(m - bound), key)) {
                        bound *= 2;

                    }
                    l = m - std::min<int>(bound, size);
                    r = m - bound / 2;
                } else {
                    int size = data_capacity_ - m;
                    while (bound < size && key_less(ALEX_DATA_NODE_KEY_AT(m + bound), key)) {
                        bound *= 2;

                    }
                    l = m + bound / 2;
                    r = m + std::min<int>(bound, size);
                }
                return binary_search_lower_bound(l, r, key);
            }

            // Searches for the first position no less than key in range [l, r)
            // https://stackoverflow.com/questions/6443569/implementation-of-c-lower-bound
            // Returns position in range [l, r]
            template <class K>
            inline int binary_search_lower_bound(int l, int r, const K& key) const {
                while (l < r) {
                    int mid = l + (r - l) / 2;
                    if (key_greaterequal(ALEX_DATA_NODE_KEY_AT(mid), key)) {
                        r = mid;
                    } else {
                        l = mid + 1;
                    }
                }
                return l;
            }


            // Insert key into pos. The caller must guarantee that pos is a gap.
            int insert_element_at(const key_type& key, key_type* payload, int pos) {
                key_slots_[pos] = key;
//                auto k1 = key_slots_[pos+1];
//            payload_slots_[pos] = payload;
                pointer_slots_[pos].second = payload;
                pointer_slots_[pos].first = true;

                set_bit(pos);

//                if (key > maxkey){
//                    auto insert_pos = pos+1;
//                    while (insert_pos < data_capacity_ && !check_exists(insert_pos)) {
//                        ALEX_DATA_NODE_KEY_AT(insert_pos) = key+1;
//                        insert_pos++;
//                    }
//                }

                // Overwrite preceding gaps until we reach the previous element
                pos--;
//            auto check_v =  check_exists(pos);   //by chao
                while (pos >= 0 && !check_exists(pos)) {
                    ALEX_DATA_NODE_KEY_AT(pos) = key;
                    pos--;
                }

//                if (pos < 0){
////                    COUT_THIS("Thank You, my Lord");
//                    return true;
//                }
                return pos;
            }

            // Insert key into pos. The caller must guarantee that pos is a gap.
            int insert_element_at(const key_type& key, bool flag, void* payload, int pos) {
                key_slots_[pos] = key;
//                auto k1 = key_slots_[pos+1];
//            payload_slots_[pos] = payload;
                pointer_slots_[pos].second = payload;
                pointer_slots_[pos].first = flag;

                set_bit(pos);

                // Overwrite preceding gaps until we reach the previous element
                pos--;
//            auto check_v =  check_exists(pos);   //by chao
                while (pos >= 0 && !check_exists(pos)) {
                    ALEX_DATA_NODE_KEY_AT(pos) = key;
                    pos--;
                }

//                if (pos < 0){
////                    COUT_THIS("Thank You, my Lord");
//                    return true;
//                }
//                return false;
                return pos;
            }

            // Insert key into pos, shifting as necessary in the range [left, right)
            // Returns the actual position of insertion
            int insert_using_shifts(const key_type& key, value_type payload, int pos) {
                // Find the closest gap
                int gap_pos = closest_gap(pos);
                set_bit(gap_pos);
                if (gap_pos >= pos) {
                    for (int i = gap_pos; i > pos; i--) {
#if ALEX_DATA_NODE_SEP_ARRAYS
                        key_slots_[i] = key_slots_[i - 1];

//                    payload_slots_[i] = payload_slots_[i - 1];
                        pointer_slots_[i].second = pointer_slots_[i - 1].second;
                        pointer_slots_[i].first = pointer_slots_[i - 1].first;

//                        if (key_slots_[i - 1] == 13944583693151){
//                            auto flag = pointer_slots_[i - 1].first;
//                            auto e = check_exists(i-1);
//                            auto e2 = check_exists(i);
//                            COUT_THIS("Jesus, please come! key_slots_[i + 1] == 17533265041071");
//                        }
                        if(check_exists(i-1)){
                            set_bit(i);
                            unset_bit(i-1);
                        }

#else
                        data_slots_[i] = data_slots_[i - 1];
#endif
                    }

                    int cflag = insert_element_at(key, payload, pos);
//                    num_shifts_ += gap_pos - pos;
//                    return pos;
                    return cflag;
                } else {
                    for (int i = gap_pos; i < pos - 1; i++) {
#if ALEX_DATA_NODE_SEP_ARRAYS

                        key_slots_[i] = key_slots_[i + 1];
//                    payload_slots_[i] = payload_slots_[i + 1];
                        pointer_slots_[i].second = pointer_slots_[i + 1].second;
                        pointer_slots_[i].first = pointer_slots_[i + 1].first;

                        if(check_exists(i+1)){
                            set_bit(i);
                            unset_bit(i+1);
                        }

#else
                        data_slots_[i] = data_slots_[i + 1];
#endif
                    }

                    int cflag = insert_element_at(key, payload, pos - 1);
//                    num_shifts_ += pos - gap_pos - 1;
//                    return pos - 1;
                    return cflag;
                }
            }

#if ALEX_USE_LZCNT
            // Returns position of closest gap to pos
        // Returns pos if pos is a gap
        int closest_gap(int pos) const {
            pos = std::min(pos, data_capacity_ - 1);
            int bitmap_pos = pos >> 6;
            int bit_pos = pos - (bitmap_pos << 6);
            // by chao, _mm_popcnt_u64 求一个整数中二进制 1 的个数
            // Count the number of bits set to 1 in unsigned 64-bit integer a, and return that count in dst.
            if (bitmap_[bitmap_pos] == static_cast<uint64_t>(-1) ||
                (bitmap_pos == bitmap_size_ - 1 &&
                 _mm_popcnt_u64(bitmap_[bitmap_pos]) ==
                 data_capacity_ - ((bitmap_size_ - 1) << 6))) {
                // no gaps in this block of 64 positions, start searching in adjacent
                // blocks
                int left_bitmap_pos = 0;
                int right_bitmap_pos = ((data_capacity_ - 1) >> 6);  // inclusive
                int max_left_bitmap_offset = bitmap_pos - left_bitmap_pos;
                int max_right_bitmap_offset = right_bitmap_pos - bitmap_pos;
                int max_bidirectional_bitmap_offset =
                        std::min<int>(max_left_bitmap_offset, max_right_bitmap_offset);
                int bitmap_distance = 1;
                while (bitmap_distance <= max_bidirectional_bitmap_offset) {
                    uint64_t left_bitmap_data = bitmap_[bitmap_pos - bitmap_distance];
                    uint64_t right_bitmap_data = bitmap_[bitmap_pos + bitmap_distance];
                    if (left_bitmap_data != static_cast<uint64_t>(-1) &&
                        right_bitmap_data != static_cast<uint64_t>(-1)) {
                        int left_gap_pos = ((bitmap_pos - bitmap_distance + 1) << 6) -
                                           static_cast<int>(_lzcnt_u64(~left_bitmap_data)) -
                                           1;
                        int right_gap_pos = ((bitmap_pos + bitmap_distance) << 6) +
                                            static_cast<int>(_tzcnt_u64(~right_bitmap_data));
                        if (pos - left_gap_pos <= right_gap_pos - pos ||
                            right_gap_pos >= data_capacity_) {
                            return left_gap_pos;
                        } else {
                            return right_gap_pos;
                        }
                    } else if (left_bitmap_data != static_cast<uint64_t>(-1)) {
                        int left_gap_pos = ((bitmap_pos - bitmap_distance + 1) << 6) -
                                           static_cast<int>(_lzcnt_u64(~left_bitmap_data)) -
                                           1;
                        // also need to check next block to the right
                        if (bit_pos > 32 && bitmap_pos + bitmap_distance + 1 < bitmap_size_ &&
                            bitmap_[bitmap_pos + bitmap_distance + 1] !=
                            static_cast<uint64_t>(-1)) {
                            int right_gap_pos =
                                    ((bitmap_pos + bitmap_distance + 1) << 6) +
                                    static_cast<int>(
                                            _tzcnt_u64(~bitmap_[bitmap_pos + bitmap_distance + 1]));
                            if (pos - left_gap_pos <= right_gap_pos - pos ||
                                right_gap_pos >= data_capacity_) {
                                return left_gap_pos;
                            } else {
                                return right_gap_pos;
                            }
                        } else {
                            return left_gap_pos;
                        }
                    } else if (right_bitmap_data != static_cast<uint64_t>(-1)) {
                        int right_gap_pos = ((bitmap_pos + bitmap_distance) << 6) +
                                            static_cast<int>(_tzcnt_u64(~right_bitmap_data));
                        if (right_gap_pos < data_capacity_) {
                            // also need to check next block to the left
                            if (bit_pos < 32 && bitmap_pos - bitmap_distance > 0 &&
                                bitmap_[bitmap_pos - bitmap_distance - 1] !=
                                static_cast<uint64_t>(-1)) {
                                int left_gap_pos =
                                        ((bitmap_pos - bitmap_distance) << 6) -
                                        static_cast<int>(
                                                _lzcnt_u64(~bitmap_[bitmap_pos - bitmap_distance - 1])) -
                                        1;
                                if (pos - left_gap_pos <= right_gap_pos - pos ||
                                    right_gap_pos >= data_capacity_) {
                                    return left_gap_pos;
                                } else {
                                    return right_gap_pos;
                                }
                            } else {
                                return right_gap_pos;
                            }
                        }
                    }
                    bitmap_distance++;
                }
                if (max_left_bitmap_offset > max_right_bitmap_offset) {
                    for (int i = bitmap_pos - bitmap_distance; i >= left_bitmap_pos; i--) {
                        if (bitmap_[i] != static_cast<uint64_t>(-1)) {
                            return ((i + 1) << 6) - static_cast<int>(_lzcnt_u64(~bitmap_[i])) -
                                   1;
                        }
                    }
                } else {
                    for (int i = bitmap_pos + bitmap_distance; i <= right_bitmap_pos; i++) {
                        if (bitmap_[i] != static_cast<uint64_t>(-1)) {
                            int right_gap_pos =
                                    (i << 6) + static_cast<int>(_tzcnt_u64(~bitmap_[i]));
                            if (right_gap_pos >= data_capacity_) {
                                return -1;
                            } else {
                                return right_gap_pos;
                            }
                        }
                    }
                }
                return -1;
            } else {
                // search within block of 64 positions
                uint64_t bitmap_data = bitmap_[bitmap_pos];
                int closest_right_gap_distance = 64;
                int closest_left_gap_distance = 64;
                // Logically gaps to the right of pos, in the bitmap these are gaps to the
                // left of pos's bit
                // This covers the case where pos is a gap
                // For example, if pos is 3, then bitmap '10101101' -> bitmap_right_gaps
                // '01010000'
                uint64_t bitmap_right_gaps = ~(bitmap_data | ((1ULL << bit_pos) - 1));
                if (bitmap_right_gaps != 0) {
                    closest_right_gap_distance =
                            static_cast<int>(_tzcnt_u64(bitmap_right_gaps)) - bit_pos;
                } else if (bitmap_pos + 1 < bitmap_size_) {
                    // look in the next block to the right
                    closest_right_gap_distance =
                            64 + static_cast<int>(_tzcnt_u64(~bitmap_[bitmap_pos + 1])) -
                            bit_pos;
                }
                // Logically gaps to the left of pos, in the bitmap these are gaps to the
                // right of pos's bit
                // For example, if pos is 3, then bitmap '10101101' -> bitmap_left_gaps
                // '00000010'
                uint64_t bitmap_left_gaps = (~bitmap_data) & ((1ULL << bit_pos) - 1);
                if (bitmap_left_gaps != 0) {
                    closest_left_gap_distance =
                            bit_pos - (63 - static_cast<int>(_lzcnt_u64(bitmap_left_gaps)));
                } else if (bitmap_pos > 0) {
                    // look in the next block to the left
                    closest_left_gap_distance =
                            bit_pos + static_cast<int>(_lzcnt_u64(~bitmap_[bitmap_pos - 1])) +
                            1;
                }

                if (closest_right_gap_distance < closest_left_gap_distance &&
                    pos + closest_right_gap_distance < data_capacity_) {
                    return pos + closest_right_gap_distance;
                } else {
                    return pos - closest_left_gap_distance;
                }
            }
        }
#else
            // A slower version of closest_gap that does not use lzcnt and tzcnt
            // Does not return pos if pos is a gap
            int closest_gap(int pos) const {
                int max_left_offset = pos;  // what's Meaning?  最多可以向左 shift 的数量
                int max_right_offset = data_capacity_ - pos - 1;  // 最多可以向右shift 的数量
                int max_bidirectional_offset =
                        std::min<int>(max_left_offset, max_right_offset);


                int distance = 1;
                while (distance <= max_bidirectional_offset) {
                      if (!check_exists(pos - distance)) {
                        return pos - distance;
                      }
                    auto insert_pos = pos + distance;
                    if (!check_exists(insert_pos)) {
                        return insert_pos;
                    }
                    distance++;
                }
                if (max_left_offset > max_right_offset) {
                    for (int i = pos - distance; i >= 0; i--) {
                        if (!check_exists(i)) return i;
                    }
                } else {
                    for (int i = pos + distance; i < data_capacity_; i++) {
                        if (!check_exists(i)) return i;
                    }
                }
                return -1;
            }
#endif


            // Insert key into pos, shifting as necessary in the range [left, right)
            // Returns the actual position of insertion
            bool insert_using_shifts(const key_type& key, bool flag, void* payload, int pos) {
                // Find the closest gap
                int gap_pos = closest_gap(pos);
                set_bit(gap_pos);
                if (gap_pos >= pos) {
                    for (int i = gap_pos; i > pos; i--) {
#if ALEX_DATA_NODE_SEP_ARRAYS
                        key_slots_[i] = key_slots_[i - 1];

//                    payload_slots_[i] = payload_slots_[i - 1];
                        pointer_slots_[i].second = pointer_slots_[i - 1].second;
                        pointer_slots_[i].first = pointer_slots_[i - 1].first;

//                        if (key_slots_[i - 1] == 13944583693151){
//                            auto flag = pointer_slots_[i - 1].first;
//                            auto e = check_exists(i-1);
//                            auto e2 = check_exists(i);
//                            COUT_THIS("Jesus, please come! key_slots_[i + 1] == 17533265041071");
//                        }
                        if(check_exists(i-1)){
                            set_bit(i);
                            unset_bit(i-1);
                        }

#else
                        data_slots_[i] = data_slots_[i - 1];
#endif
                    }

                    bool cflag = insert_element_at(key,flag, payload, pos);
                    num_shifts_ += gap_pos - pos;
//                    return pos;
                    return cflag;
                } else {
                    for (int i = gap_pos; i < pos - 1; i++) {
#if ALEX_DATA_NODE_SEP_ARRAYS

                        key_slots_[i] = key_slots_[i + 1];
//                    payload_slots_[i] = payload_slots_[i + 1];
                        pointer_slots_[i].second = pointer_slots_[i + 1].second;
                        pointer_slots_[i].first = pointer_slots_[i + 1].first;

                        if(check_exists(i+1)){
                            set_bit(i);
                            unset_bit(i+1);
                        }

#else
                        data_slots_[i] = data_slots_[i + 1];
#endif
                    }

                    bool cflag = insert_element_at(key, flag, payload, pos - 1);
                    num_shifts_ += pos - gap_pos - 1;
//                    return pos - 1;
                    return cflag;
                }
            }



            // Model-based inserts
            forceinline int  model_insert(key_type key, key_type* payload, int predict_pos) {
                // predict_pos is the position that use the model to predict the position of key to be inserted
//                predict_pos = std::max<int>(std::min<int>(predict_pos, data_capacity_ - 1), 0);
                predict_pos = std::min<int>(predict_pos, data_capacity_ - 1);
                // Insert
//                bool cflag = false;
                ++num_keys_;
                if (num_keys_){
//                    if (!check_exists(predict_pos) && !check_exists(std::max<int>(0,predict_pos-1)))
//                        cflag = insert_element_at(key, payload, predict_pos);
//                    else if (!check_exists(predict_pos) && key_less(ALEX_DATA_NODE_KEY_AT(predict_pos-1), key) ){
//                        cflag = insert_element_at(key, payload, predict_pos);
//                    }
                    if (!check_exists(predict_pos) && key_lessequal(ALEX_DATA_NODE_KEY_AT(std::max<int>(0,predict_pos-1)), key) ){

                        return insert_element_at(key, payload, predict_pos);
                    }
                    else{
                        std::pair<int, int> positions = find_insert_position(key,predict_pos);
//                        int upper_bound_pos = positions.second;
//                        if (!allow_duplicates && upper_bound_pos > 0 &&
//                            key_equal(ALEX_DATA_NODE_KEY_AT(upper_bound_pos - 1), key)) {
//                            return {-1, upper_bound_pos - 1};
//                        }
                        int insertion_position = positions.first;
                        if (insertion_position < data_capacity_ &&
                            !check_exists(insertion_position)) {
                            bool cflag = insert_element_at(key, payload, insertion_position);
                            return cflag;
                        } else {
                            return insert_using_shifts(key, payload, insertion_position);
                        }
                    }
                }
                else{
//                    cflag = insert_element_at(key, payload, predict_pos);
                    return insert_element_at(key, payload, predict_pos);
                }

                /*
                // 什么时候会发生 lru node change
                if (++num_keys_ == data_capacity_){
//                    COUT_THIS("Jesus, please come! retrain this leaf, num_keys ++ == data_capacity_ ");
                    return {true,cflag};
                }
                else{
                    return {false,cflag};
                }
                */
//                return cflag;

            }

            // Model-based inserts, during retrain a single
            forceinline void model_insert(key_type key, bool flag, void* payload, int predict_pos) {
                // predict_pos is the position that use the model to predict the position of key to be inserted
                predict_pos = std::max<int>(std::min<int>(predict_pos, data_capacity_ - 1), 0);
                // Insert
//                bool cflag = false;
                if (num_keys_){
                    if (!check_exists(predict_pos) && !check_exists(std::max<int>(0,predict_pos-1)))
                        insert_element_at(key,flag, payload, predict_pos);
                    else if (!check_exists(predict_pos) && key_less(ALEX_DATA_NODE_KEY_AT(predict_pos-1), key) ){
                        insert_element_at(key, flag, payload, predict_pos);
                    }
                    else{
                        std::pair<int, int> positions = find_insert_position(key,predict_pos);
//                        int upper_bound_pos = positions.second;
//                        if (!allow_duplicates && upper_bound_pos > 0 &&
//                            key_equal(ALEX_DATA_NODE_KEY_AT(upper_bound_pos - 1), key)) {
//                            return {-1, upper_bound_pos - 1};
//                        }
                        int insertion_position = positions.first;
                        if (insertion_position < data_capacity_ &&
                            !check_exists(insertion_position)) {
                            insert_element_at(key, flag, payload, insertion_position);
                        } else {
                            insertion_position =
                                    insert_using_shifts(key, flag, payload, insertion_position);
                        }
                    }
                }
                else{
                    insert_element_at(key, flag, payload, predict_pos);
                }
                // 什么时候会发生 lru node change
                ++num_keys_;
//                if (++num_keys_ == data_capacity_){
////                    COUT_THIS("Jesus, please come! retrain this leaf, num_keys ++ == data_capacity_ ");
//                    return {true,cflag};
//                }
//                else{
//                    return {false,cflag};
//                }

            }



            // inner level rebuild 时，unlink leaf piece 时插入buffer key, 因为buffer key 已经是有序的，所以直接向后插入
            forceinline pair<bool,bool>  insert(key_type key, bool flag, void* payload) {
                this->slotkey.emplace_back(key);

                this->slotdata.emplace_back( payload);
                this->locbitmap.emplace_back(flag);
                ++this->intercept;
            }


            // 判断key 是否在当前buffer 中
            inline pair<bool,unsigned long> is_in_leaf(const key_type &key,const int error){
                pair<bool,unsigned long> is_in(false,0);
                if (key<=this->endkey ){
                    auto pos = (*this)(key);

                    int slotlo = PGM_SUB_EPS(pos, error + 1) ;  //,(--itlevel)->size()
                    int slothi = PGM_ADD_EPS(pos, error,(this->slotkey.size()-1));
                    for (; slotlo<=slothi;slotlo++){

                        if (key == this->slotkey[slotlo])
                        {   is_in.first = true;
                            is_in.second = slotlo;
                            break;}
                    }
//            cout<< " happy new year! Jesus~~~~" << endl;
                    return is_in;
                }
                else
                    return is_in;
            }

            // 找到key 在当前buffer 中，最后一个小于或等于key 的position
            inline unsigned long find_upper_in_leaf(const key_type &key, int predict_pos){
//                pair<bool,unsigned long> is_in(false,0);
                predict_pos = std::max<int>(std::min<int>(predict_pos, data_capacity_ - 1), 0);
//                predict_pos = std::max<int>(std::min<int>(predict_pos, data_capacity_ - 1), 0);
                int pos = exponential_search_upper_bound(predict_pos, key);

//                 is_in.first = true;
//                 is_in.second = pos;

                 return pos;
            }

        };

        struct Innerpiece{

            key_type startkey;
            double slope;
            double intercept;

            inline Innerpiece() = default;

            inline Innerpiece(key_type startkey,double slope, double intercept):
                    startkey(startkey),slope(slope),intercept(intercept){};

            inline explicit Innerpiece(double slope, double intercept):
                    startkey(),slope(slope),intercept(intercept){};


            inline explicit Innerpiece(const typename internal::insertPWLF<key_type, int>::CanonicalSegment &cs):
                    startkey(cs.get_first_x()){
                auto[cs_slope, cs_intercept] = cs.get_floating_point_segment(startkey);
                if (cs_intercept > std::numeric_limits<decltype(intercept)>::max())
                    throw std::overflow_error("Change the type of Segment::intercept to int64");
                slope = cs_slope;
                intercept = cs_intercept;
            };

            friend inline bool operator<(const Innerpiece &p, const key_type &key) { return p.startkey < key; }
            friend inline bool operator<(const key_type &key, const Innerpiece &p) { return key < p.startkey; }
            friend inline bool operator<(const Innerpiece &p1, const Innerpiece &p2) { return p1.startkey < p2.startkey; }

            operator key_type() { return startkey; };

            /**
         * Returns the approximate position of the specified key.
         * @param k the key whose position must be approximated
         * @return the approximate position of the specified key
         */
            forceinline size_t operator()(const key_type &key) const {
                auto pos = int64_t(slope * (key - startkey)) + intercept;
                return pos > 0 ? size_t(pos) : 0ull;
            }

            inline void update(typename internal::insertPWLF<key_type, int>::CanonicalSegment &cs){
                startkey = cs.first;
                auto[cs_slope, cs_intercept] = cs.get_floating_point_segment(startkey);
                if (cs_intercept > std::numeric_limits<decltype(intercept)>::max())
                    throw std::overflow_error("Change the type of Segment::intercept to int64");
                slope = cs_slope;
                intercept = cs_intercept;
            }


        };

        // 声明 innerlevel 结构体，每个innerllevel 包含一个 vector<Innerpiece>,一个 opt， 一个 pos；
        struct Innerlevel{
            std::vector<Innerpiece* > innerpieces;
//        vector<internal::insertPWLF<key_type, int>*> opt;
            internal::insertPWLF<key_type, int>* opt;
            unsigned int pos = 0;
            unsigned int nextpos = 0;
            Innerpiece* i_innerpiece = new Innerpiece();  // 该innerlevel 的最后一个inner piece
        };


        struct Leaflevel{
            std::vector<Leafpiece* > leafpieces;
//        vector<internal::insertPWLF<key_type, int>*> opt;
            internal::insertPWLF<key_type, int>* opt;
            unsigned int pos = 0;
        };

        unsigned int totalnum;
        double buffer_ratio;
        unsigned int retrain = 0;  // 统计model retrain 的次数, 当该leaf 的buffer piece 达到阈值时
        int newmodels = 0;  // 由于retrain带来的多出来的new models的个数。只有当retrain 时models 个数增加，才会导致leaves 个数变化
        unsigned rebuild_inner = 0; // 记录inner level trebuild 的次数
        unsigned long int inkeynum=0;
        unsigned long int exkeynum=0;
        int maxbsize = 0;
        int minbsize = 9999999;  // 2147483647
        int valuesize = 0;
        int Error ;
        double fill_factor;
        unsigned long key_in_buff = 0;
        unsigned long reserve_buff = 0;
        int ErrorRecursive ;
        int leafsplit = 0;
        Leafpiece *m_transleaf = NULL;
        Leafpiece *m_tailleaf = NULL;   // tailpiece, 用于append increasing keys
//        Innerpiece *i_innerpiece = NULL;
        queue< pair<key_type,unsigned int> > buffqueue;
        film_stats vstats;

        Innerpiece* root;
//        typedef std::pair< std::vector<Leafpiece*>,vector<internal::insertPWLF<key_type, int>*> > leafleveltype;
        Leaflevel leaflevel;   //  the leafpieces composing the leaflevel
//        typedef std::pair< std::vector<Innerpiece>, vector<internal::insertPWLF<key_type, int>*> > innerleveltype;
        vector<Innerlevel*> innerlevels;



        /**
         * constructs an empty film
         */
//        FILMinsert() = default;

        /**
         * consttructs index on the given keys vector
         * @param keys the vector of keys to be indexed (sorted). @param payload the value of each key
         */

        FILMinsert(unsigned int n ,int err1,int err2)
        {
            totalnum =n,Error = err1, ErrorRecursive =err2;

//            internal::insertPWLF<key_type, int> opt(err1);
//            leaflevel.second.push_back(&opt) ;
        }

        ~FILMinsert() {
            // delete leaf piece
        }
        static void delete_leaf(Leafpiece* node) {
            if (node == NULL) {
                return;
            } else{
//                vector<key_type>().swap(node->slotkey);
//                vector<bool>().swap(node->locbitmap);
//                vector<void*>().swap(node->slotdata);
                // release intrachain
                while (node->buffer_piece != NULL){
                    node = node->buffer_piece;
                    if (node->buffer_piece == NULL){
                        sortPiece* buffer = (sortPiece*) node;
                        delete buffer;
                        break;}
                    vector<key_type>().swap(node->slotkey);
                    vector<bool>().swap(node->locbitmap);
//                    vector<void*>().swap(node->slotdata);
                    node->intrachain.deletelru();
                }
                node->buffer_piece = NULL;
                delete node->buffer_piece;
//                node->intrachain.deletelru();
//                for (int i = 0; i < node->intrachain.size; i ++){
//                    delete node->intrachain[i];
//                }
//                adalru::localLRU <lruOff_type,key_type* >  intrachain;    //unsigned short int
//                delete node;
                node = NULL;
//                data_node_allocator().deallocate(static_cast<data_node_type*>(node), 1);
            }
        }

        void delete_inner(Innerpiece* node) {

            if (node == NULL)
                return;
            else
                delete[] node;
            node = nullptr;
//                data_node_allocator().deallocate(static_cast<data_node_type*>(node), 1);
        }

        void release(){
            for ( unsigned int leafi = 0; leafi < leaflevel.leafpieces.size();leafi++){
                auto leaf_it = leaflevel.leafpieces[leafi];
                delete_leaf(leaf_it);
            }
            vector<Leafpiece*>().swap(leaflevel.leafpieces);
            delete leaflevel.opt;
            for (unsigned int leveli = 0; leveli <innerlevels.size();leveli ++ ){
                auto level_it = innerlevels[leveli];
                innerlevels[leveli]->i_innerpiece = NULL;
                for (unsigned int inneri = 0; inneri < level_it->innerpieces.size();inneri ++ ){
                    auto inner_it = level_it->innerpieces[inneri];
                    delete_inner(inner_it);
                }
//                level_it->innerpieces.clear();
                vector<Innerpiece*>().swap(level_it->innerpieces);
                delete level_it->opt;
            }
            innerlevels.clear();

            inkeynum = 0;
            exkeynum = 0;
            m_transleaf = NULL;
            m_tailleaf = NULL;
        }

        void release_inner(){

            for (unsigned int leveli = 0; leveli <innerlevels.size();leveli ++ ){
                auto level_it = innerlevels[leveli];
                innerlevels[leveli]->i_innerpiece = NULL;
                for (unsigned int inneri = 0; inneri < level_it->innerpieces.size();inneri ++ ){
                    auto inner_it = level_it->innerpieces[inneri];
                    delete_inner(inner_it);
                }
//                level_it->innerpieces.clear();
                vector<Innerpiece*>().swap(level_it->innerpieces);
                delete level_it->opt;
            }
            innerlevels.clear();

        }

        double compute_fill_factor(){
            key_in_buff = 0;
            reserve_buff = 0;
            for (int i = 0; i < leaflevel.leafpieces.size();i++)
            {
                Leafpiece * curleaf = leaflevel.leafpieces[i];
                while (curleaf->buffer_piece!= NULL)
                    curleaf = curleaf->buffer_piece;
                sortPiece* buffer = (sortPiece*) curleaf;
                key_in_buff += buffer->num_keys_;
                reserve_buff += buffer->data_capacity_;
            }
            fill_factor = (double)key_in_buff/reserve_buff;
        }

        /*
         * the method to construct film, that build FILM
         */
        template<  typename lru_type>
        inline void internal_level_rebuild (unsigned short int error_recursize,lru_type interchain){
            // 释放原来的inner levels
            this->release_inner();
            Innerlevel *innerlevel0 = new Innerlevel;
            std::pair<size_t, vector<key_type> > n_parts = internal::make_segmentation(this, error_recursize, &leaflevel, innerlevel0, leaflevel.leafpieces[0]->startkey,interchain);
            innerlevels.emplace_back(innerlevel0);
            while (n_parts.first >1 ){
                Innerlevel *innerlevel = new Innerlevel;
                n_parts = internal::make_segmentation(this, error_recursize,n_parts.second,innerlevel);   //<key_type,Leafpiece>
                //                std::cout<< "You are my refuge, from ever to ever"<<endl;

                innerlevels.emplace_back(innerlevel);
            }
            root = innerlevels.back()->innerpieces[0];
        }

        inline std::vector<key_type>  append_one(size_t error,std::vector<key_type> keys,unsigned int k){

            std::vector<key_type> startkeys;

            for (size_t i = 0; i < keys.size(); ++i) {
                pair<key_type,unsigned int> p(keys[i],innerlevels[k]->nextpos++) ;   // i 为 pos
                ++(innerlevels[k]->pos);
                if (!innerlevels[k]->opt->add_point(p.first, p.second)) {  // 如果inner level  不满足error 了，那么再创建一个innerpiece
                    // 当前 innerpiece 不再满足，需要创建new inner piece 并判断该 innerlevel 的上一层level 是否需要更新
                    auto a = innerlevels[k]->opt->get_segment();
//                    Innerpiece* innerpiece;//
                    innerlevels[k]->i_innerpiece->update(a);
                    if (innerlevels[k]->pos > 2)
                        innerlevels[k]->innerpieces.pop_back();
                    innerlevels[k]->innerpieces.emplace_back(innerlevels[k]->i_innerpiece);
//                    auto newinner = innerlevels[k]->innerpieces.back();
//                    cout<< "i need You, my lovely Lord, why there is nan" << endl;

                    // 首先在该层创建一个 new innerpiece， 更新该innerpiece，再递归向向上
                    Innerpiece *innerpiece = new Innerpiece();
                    innerlevels[k]->i_innerpiece = innerpiece;
                    innerlevels[k]->pos = 0;
                    innerlevels[k]->nextpos -= 2;
                    delete innerlevels[k]->opt;
//                    innerlevels[k]->opt.pop_back();
                    internal::insertPWLF<key_type, int> *inneropt = new internal::insertPWLF<key_type, int>(error);
                    innerlevels[k]->opt = inneropt;
                    if (k==0){
//                        auto aaaaa = leaflevel.leafpieces.size()-2;
                        startkeys.emplace_back(leaflevel.leafpieces[leaflevel.leafpieces.size()-1]->startkey);
                    }

                    else{
                        startkeys.emplace_back(innerlevels[k-1]->innerpieces[innerlevels[k-1]->innerpieces.size()-2]->startkey);
                    }
                    startkeys.emplace_back(p.first);
                    auto rr = append_one(error,startkeys,k);

                    if (innerlevels.back()->innerpieces.size() > 1)
                    {
                        startkeys.clear();
                        startkeys.emplace_back(a.first);
                        startkeys.emplace_back(p.first);
//                        Innerpiece innerpiece;// 创建parent piece
                        internal::insertPWLF<key_type, int> *inneropt = new internal::insertPWLF<key_type, int>(error);
//                    std::pair< std::vector<typename filmtype:: Innerpiece>, std::vector<internal::insertPWLF<key_type, int>*> > *innerlevel =
//                            new std::pair< std::vector<typename filmtype:: Innerpiece>, std::vector<internal::insertPWLF<key_type, int>*> >;
                        Innerlevel *innerlevel = new Innerlevel;
                        innerlevels.emplace_back(innerlevel);
                        innerlevels.back()->opt = inneropt ;
                        auto rr = append_one(error,startkeys,k+1);
//                    cout<< "Jesus, i need You !"  << endl;

                        return startkeys ;
                    }
                    else if (innerlevels.size() > 1 && (k != innerlevels.size()-1) )   // 如上为 由于创建了new innner piece，导致了new innerlevel，如下为，虽然创建了new innerpiece，但只需要更新上层的innner level
                    {
                        // 更新上层的最后一个inner piece
                        startkeys.pop_back();
                        auto rr = append_one(error,startkeys,k+1);
                        startkeys.clear();
//                    cout << "thank You, my Lord! i need You!"<<endl;
//                        startkeys.emplace_back(a.first);
                        return startkeys ;
                    }
                    else if (innerlevels.size() > 1 && (k == innerlevels.size()-1) )   // 如上为 由于创建了new innner piece，导致了new innerlevel，如下为，虽然创建了new innerpiece，但只需要更新上层的innner level
                    {
//                        cout << "thank You, my Lord! i need You!"<<endl;
                        startkeys.clear();
                        return startkeys ;
                    }
                    else{
                        startkeys.clear();
                    }
                    a = innerlevels[k]->opt->get_segment();
//                    if (k == 0 && innerlevels[0]->innerpieces.size() == 6){
//                        cout << "Jesus, please come!" <<endl;
//                    }
                    innerlevels[k]->i_innerpiece->update(a);
                    if (innerlevels[k]->pos > 2)
                        innerlevels[k]->innerpieces.pop_back();
                    innerlevels[k]->innerpieces.emplace_back(innerlevels[k]->i_innerpiece);
//                    if (k == 0 && innerlevels[0]->innerpieces.size() == 6){
//                        cout << "Jesus, please come!" <<endl;
//                    }

                    startkeys.emplace_back(a.first);
                    return startkeys;
                }
            }

            auto a = innerlevels[k]->opt->get_segment();

            innerlevels[k]->i_innerpiece->update(a);
            if (innerlevels[k]->pos > 2)
                innerlevels[k]->innerpieces.pop_back();
            innerlevels[k]->innerpieces.emplace_back(innerlevels[k]->i_innerpiece);

//            if (k == 0 && innerlevels[0]->innerpieces.size() == 6){
//                cout << "Jesus, please come!" <<endl;
//            }

            startkeys.emplace_back(a.first);
            return startkeys ;

        }



        template<class lru_type>   // insert one key on an existing index
        forceinline void append_one(key_type key,key_type* payload,unsigned int error,lru_type interchain){
            std::vector<key_type> startkeys;
            inkeynum++;
//            if (key == 60721668)
//                cout << "i need You, my Lord! " << endl;
            pair<key_type,unsigned int> p(key,leaflevel.opt->points_in_hull) ;

            if (leaflevel.opt->append_point(p.first, p.second)) {
                m_tailleaf->slotkey.emplace_back(p.first);
                m_tailleaf->endkey = p.first;
                m_tailleaf->slotdata.emplace_back( m_tailleaf->intrachain.put(p.second,payload));
                m_tailleaf->locbitmap.emplace_back( true);
            }
            else
            {
                auto a = leaflevel.opt->get_segment(m_tailleaf->endkey); // 将生成的new leaf piece插入到leaflevel 中
                m_tailleaf->update(a);
                interchain->put(m_tailleaf->startkey,m_tailleaf);
                auto buffersize = ceil(m_tailleaf->slotkey.size()*this->buffer_ratio);   //0.429, 0.25, 0.667
                if (buffersize < minbsize)
                    minbsize = buffersize;
                else if (buffersize > maxbsize)
                    maxbsize = buffersize;
                m_tailleaf->buffer_piece->slope = buffersize;
                m_tailleaf->buffer_piece->slotkey.reserve(buffersize);
                m_tailleaf->buffer_piece->slotdata.reserve(buffersize);
                m_tailleaf->buffer_piece->locbitmap.reserve(buffersize);
                // 这里是初始化 parent piece 的first key 和 last key
                if (innerlevels.size() == 0){
                    startkeys.emplace_back( m_tailleaf->startkey);
                    startkeys.emplace_back( p.first);
                    Innerpiece innerpiece;// 创建parent piece
                    internal::insertPWLF<key_type, int> *inneropt = new internal::insertPWLF<key_type, int>(error);
//                    std::pair< std::vector<typename filmtype:: Innerpiece>, std::vector<internal::insertPWLF<key_type, int>*> > *innerlevel =
//                            new std::pair< std::vector<typename filmtype:: Innerpiece>, std::vector<internal::insertPWLF<key_type, int>*> >;
                    Innerlevel *innerlevel = new Innerlevel;
                    innerlevels.emplace_back(innerlevel);
                    innerlevels[0]->opt = inneropt ;
//                    cout << " my Lord, Jesus, please have pity on me"<< endl;
                }
                else{
                    // 从innerlevel 的最底层到root 层，判断是否需要更新
                    startkeys.emplace_back( p.first);  //p.first is the break_key, and also the startkey of the next leaf piece
//                    cout << "my lovely Lord, i trust in You!" << endl;
                }
                // use the
                auto rr = append_one(error,startkeys,0);
                startkeys.clear();
//                cout << "Jesus, i need You!!"<< endl;

                m_tailleaf = new Leafpiece;
                internal::insertPWLF<key_type, int> *leafopt = new internal::insertPWLF<key_type, int> (error);
                delete leaflevel.opt;
                leaflevel.opt = leafopt ;
                buffersize = (maxbsize+minbsize)/2;
                m_tailleaf->buffer_piece = new sortPiece(buffersize);
//                buffersize = ceil(8192*2*0.667);  //0.429
                m_tailleaf->buffer_piece->slope = buffersize;
                m_tailleaf->buffer_piece->slotkey.reserve(buffersize);
                m_tailleaf->buffer_piece->slotdata.reserve(buffersize);
                m_tailleaf->buffer_piece->locbitmap.reserve(buffersize);
                m_tailleaf->slotkey.reserve(8192*2);
                leaflevel.opt->append_point(p.first, 0);
                m_tailleaf->slotdata.emplace_back( m_tailleaf->intrachain.put(0,payload));
                m_tailleaf->slotkey.push_back(p.first);
                m_tailleaf->locbitmap.emplace_back( true);
                m_tailleaf->startkey = p.first;
                m_tailleaf->endkey = p.first;
                interchain->put(m_tailleaf->startkey,m_tailleaf);

                leaflevel.leafpieces.emplace_back( m_tailleaf);
//                auto newleaf = leaflevel.leafpieces.back();
//                cout<< "i need You, my lovely Lord!" << endl;

            }
        }




        template<class lru_type>  // insert_bulk_load_one_by_one (also one by one), in the beginning, the index must be empty? maybe not  empty is okay
        inline void update_append(std::vector<key_type> keys, key_type* payload,unsigned int error,unsigned int error_recursize,lru_type &interchain ){
            // build leaf level
            // size_t error,std::vector<key_type> keys,std::vector<key_type> payload,
            //filmtype *filmada, unsigned int &inkeynum,leaf_type* m_taillea
            std::pair<size_t,std::vector<key_type>>  n_parts = internal::append_segmentation(error,keys,payload,this,inkeynum,this->m_tailleaf,interchain);   //<key_type,Leafpiece>
//            std::cout<< "You are my refuge, from ever to ever"<<endl;
            // on the basis of leaf level, build the upper level until just has one piece, that the make_segmentation just return 1.
            if (innerlevels.size()>=1)
                root = innerlevels.back()->innerpieces[0];
            this->verify_film(&vstats);
            cout<<"index append test finished,Jesus, my Lord, i need You, praise to You for ever and ever!" << endl;
        }


        //  与 range query 相关的返回结果
        struct range_res_find
        {
            /// find result flags
            bool  flags;

            /// The key to be found
            key_type    firstkey;
            key_type    lastkey;

            /// The leaf nodes that overlap with the queried range
            vector<Leafpiece *> findleafs;

            /// the location of key in leaf node
            unsigned int firstslot;
            unsigned int lastslot;

            /// Constructor of a result with a specific flag, this can also be used
            /// as for implicit conversion.

            range_res_find(bool f)
                    : flags(f), firstkey(),lastkey(), findleafs(),firstslot(),lastslot()
            { }

            range_res_find(bool f,key_type firstkey,key_type lastkey)
                    : flags(f), firstkey(firstkey), lastkey(), findleafs(),firstslot(),lastslot()
            { }


            /// Constructor with a lastkey value.
            forceinline range_res_find(bool f ,key_type firstkey,int first_loc,key_type lastkey)
                    : flags(f), firstkey(firstkey), findleafs(),firstslot(first_loc),lastkey(lastkey),lastslot()
            { }

            inline bool in_or_out(){
                return true;
            }

        };

        //  与 range query 相关的返回结果
        struct htaprange_res_find
        {

            /// The key to be found
            key_type    firstkey;
            key_type    lastkey;

            /// The leaf nodes that overlap with the queried range
            vector<Leafpiece *> findleafs;
            vector<sortPiece *> findbuffers;

            /// the location of key in leaf node
            vector<unsigned int> leafslots;
            int slotleaf0;
            int slotleaf1;

            int slotbuffer0;
            int slotbuffer1;

            /// Constructor of a result with a specific flag, this can also be used
            /// as for implicit conversion.

//            range_res_find(bool f)
//                    : flags(f), firstkey(),lastkey(), findleafs(),firstslot(),lastslot()
//            { }

            htaprange_res_find(key_type firstkey)
                    : firstkey(),lastkey(), findleafs(),findbuffers(),leafslots(),slotleaf0(),slotleaf1(),slotbuffer0(),slotbuffer1()
            { }

            htaprange_res_find(key_type firstkey,key_type lastkey,int slotleaf0,int slotbuffer0)
                    : firstkey(firstkey),lastkey(lastkey), findleafs(),findbuffers(),leafslots(),slotleaf0(slotleaf0),slotleaf1(),slotbuffer0(slotbuffer0),slotbuffer1()
            { }

            htaprange_res_find(key_type firstkey,int slotleaf0,int slotbuffer0)
                    : firstkey(firstkey),lastkey(), findleafs(),findbuffers(),leafslots(),slotleaf0(slotleaf0),slotleaf1(),slotbuffer0(slotbuffer0),slotbuffer1()
            { }

            inline bool in_or_out(){
                return true;
            }

        };

        /// film recursive find the key has much information which is needs to be return.
        struct result_find
        {
            /// find result flags
            bool find;   // 在sorted list 还是在 regular data
            // find 表示 key 是在leaf 中，还是在buffer 中； true 表示 leaf，false 表示buffer
            /// this flag indicate whether the queried data is in memory or on disk
            bool  flags;

            /// The leaf node the findkey belong to
            Leafpiece *findleaf;

            /// the location of key in leaf node
            int slot;

            int leafslot;

            /// Constructor of a result with a specific flag, this can also be used
            /// as for implicit conversion.
            inline result_find(bool f = true)
                    : find (true),flags(f), findleaf(),slot()
            { }

            /// Constructor with a lastkey value.
            inline result_find(bool find_f,bool in_or_exf ,  Leafpiece * find_leaf,int loc_slot)
                    : find(find_f), flags(in_or_exf),  findleaf(find_leaf),slot(loc_slot)
            { }

        };

        struct range_key_result_find
        {

            /// The leaf node the findkey belong to
            vector<Leafpiece *>findleafs;   // the leaf it belong to, the buffer that overlap with search range

            /// the location of key in leaf node
            vector<int> slots;   // the slot in leaf, the slot in buffer

            int leafslot;

            /// Constructor of a result with a specific flag, this can also be used
            /// as for implicit conversion.
            inline range_key_result_find(bool f = true)
                    :  findleafs(),leafslot()
            { }

            /// Constructor with a lastkey value.
            inline range_key_result_find(int loc_slot)
                    : leafslot(loc_slot)
            { }


        };

        // Searches for the first position greater than key, starting from position m
// Returns position in range [0, data_capacity]


        /*** Comparison ***/

        template<class pieces_type>
        inline int exponential_search_upper_bound(const pieces_type childpieces,unsigned long m, const key_type key) {
            // Continue doubling the bound until it contains the upper bound. Then use
            // binary search.
            int bound = 1;
            auto childsize = childpieces.size();
            if (key_lessequal(childpieces.size(),m))
                m = childsize -1;
            int l, r;  // will do binary search in range [l, r)
            auto lo =  childpieces.begin();  //,(--itlevel)->size()
//            auto hi = childpieces.begin();
            auto ckey = (*(lo+m))->startkey;   // (*leaflo)

            // 读取位置m 的值。internal node gets the start key, leaf nose gets the key stored in m
            if (key_greater((*(lo+m))->startkey, key)) {   // if the key in predicted position is greater than the query key
                int size = m;
                while (bound < size &&
                       key_greater( (*(lo+m-bound))->startkey, key)) {
//                    key_type mkey = (*(lo+m-bound))->startkey  ;// chao
//        key_greater(mkey, key);
//                    auto mkeydiff = mkey - key;  //chao
//                    auto mkey2 = key_slots_[m-bound];   //chao
                    bound *= 2;
//                    num_exp_search_iterations_++;
                }
                l = m - std::min<int>(bound, size);
                r = m - bound / 2;
            } else {    // the key in the predicted position is smaller than the query key, by chao
                int size = childpieces.size() - m;   // search in the right of the position m, by chao
                while (bound < size &&
                       key_lessequal((*(lo + m + bound))->startkey, key)) {
                    bound *= 2;
//                    num_exp_search_iterations_++;
                }
                l = m + bound / 2;
                r = m + std::min<int>(bound, size);
            }


            return binary_search_upper_bound(childpieces, l, r, key);
        }



        // Searches for the first position greater than key in range [l, r)
// https://stackoverflow.com/questions/6443569/implementation-of-c-lower-bound
// Returns position in range [l, r]
//        template <class piece_type>
//        inline int binary_search_upper_bound(piece_type &pieces,int l, int r, const key_type key) const {
//            while (l < r) {
//                int mid = l + (r - l) / 2;
//                if (key_lessequal(pieces[mid].startkey, key)) {
//                    l = mid + 1;
//                } else {
//                    r = mid;
//                }
//            }
//            return l;
//        }

        template <class piece_type>
        inline int binary_search_upper_bound(piece_type &pieces,int l, int r, const key_type key) const {
            while (l < r) {
                int mid = l + (r - l) / 2;
                auto psize = pieces.size();
                auto piece = *pieces[mid];// by chao
                if (key_lessequal((*pieces[mid]).startkey, key)) {
                    l = mid + 1;
                } else {
                    r = mid;
                }
            }
            return l;
        }

        /** search method
         * the methods used in query, @param key is the queried key, return the leafpiece responsible for the given key
         */



        forceinline result_find search_one(const key_type key) {
            // 从 root 开始，即，从root开始的结尾开始
//            struct timeval ct1, ct2;
//            double ctimeuse;
//            if (vstats.innernum > 1) {
                auto itlevel = innerlevels.end() -1; // root level 的 iterator
                // internal level use exponential search, the parameters used in exponential_search are: predict_pos, search key
                /*
                 * 1. get the predicted pos
                 * 2. get the [l,r] used in binary search or, get the [l,r]
                 * 3. use the binary search in [l,r] to get the final position
                 */
                auto predicted_pos = (*root)(key);
                auto childpieces = (*(--itlevel))->innerpieces;
//                auto search_pos = std::max(exponential_search_upper_bound(childpieces,predicted_pos, key) -1,0);
            auto search_pos =exponential_search_upper_bound(childpieces,predicted_pos, key) -1;
//            if (search_pos <0)
//                search_pos = 0;
                // 判断 找到的inner node 是否满足条件
                auto curit = (*(itlevel))->innerpieces.begin() + search_pos;
                Innerpiece* b = *curit;
//                if (b->startkey > key )
//                    cout << "i need You, my lovely Lord, b->startkey > key " << endl;
                for (; --itlevel >= innerlevels.begin();){
                    predicted_pos = (*b)(key);
                    childpieces = (*(itlevel))->innerpieces;
                    search_pos = exponential_search_upper_bound(childpieces,predicted_pos, key) -1;
//                if (search_pos <0)
//                    search_pos = 0;
                    // 判断 找到的inner node 是否满足条件
                    curit = (*(itlevel))->innerpieces.begin() + search_pos;
                }
                // find the leaf level

                b = *curit;
                predicted_pos = (*b)(key);

                search_pos = exponential_search_upper_bound(leaflevel.leafpieces,predicted_pos, key) -1;
//            if (search_pos <0)
//                search_pos = 0;
                // 预测位置和最终位置的距离
                auto curleaf = leaflevel.leafpieces.begin() + search_pos;

                Leafpiece* a = *curleaf;
                predicted_pos = (*a)(key);
                // 确定 该leaf 的哪一个链 leaf piece
                while ( a->buffer_piece->buffer_piece && key > a->endkey) {
                    a = a->buffer_piece;
                }
                if (a->buffer_piece == NULL) {
//                    gettimeofday(&ct1, NULL);
                    sortPiece *buffer = (sortPiece *) a;
                    auto resslot = predicted_pos * this->buffer_ratio;
                    auto model_res = buffer->find_key(key, resslot);

//                if (buffer->ALEX_DATA_NODE_KEY_AT(model_res) == key) {
                    result_find index_res = result_find(false, buffer->ALEX_DATA_NODE_FLAG_AT(model_res), buffer,
                                                        model_res);   // sort_list data
                    return index_res;
//                } else {
//                    cout << " i need You, my Lord, please help me! not find the key " << key << endl;
//                }
                    /*
                    auto sortpos = (*buffer)(key);
                    if (buffer->slotkey[sortpos] == key) {
                        result_find index_res = result_find(false, buffer->locbitmap[sortpos], buffer,
                                                            sortpos);   // sort_list data
                        return index_res;
                    } else {
                        cout << " i need You, my Lord, please help me! not find the key " << key << endl;
                    }
                    */
//                    gettimeofday(&ct2, NULL);
//                    ctimeuse = (ct2.tv_sec - ct1.tv_sec) + (double) (ct2.tv_usec - ct1.tv_usec) / 1000000.0;
//                    query_stats->sort_piecetime += ctimeuse;
                }
                else {
                    predicted_pos = (*a)(key);

                    int slotlo = PGM_SUB_EPS(predicted_pos, ErrorRecursive + 1);
                    slotlo = std::min<int>(slotlo,(a->slotkey.size() - 1));
                    int slothi = PGM_ADD_EPS(predicted_pos, ErrorRecursive, (a->slotkey.size() - 1));
                    int resslot = slotlo;
                    for (; resslot < slothi; resslot++) {
                        if (key <= a->slotkey[resslot]) {
                            break;
                        }
                    }
//                    gettimeofday(&ct1, NULL);
                    if (key != a->slotkey[resslot]) {
                        // find in sort_list
                        while (a->buffer_piece->buffer_piece != NULL) {
                            a = a->buffer_piece;
                        }
                        predicted_pos = (*(Leafpiece*)(*curleaf))(key);
                        sortPiece *buffer = (sortPiece *) a->buffer_piece;
                        resslot = predicted_pos * this->buffer_ratio;
                        auto model_res = buffer->find_key(key, resslot);
//                    if (buffer->ALEX_DATA_NODE_KEY_AT(model_res) == key) {
                        result_find index_res = result_find(false, buffer->ALEX_DATA_NODE_FLAG_AT(model_res), buffer,
                                                            model_res);   // sort_list data
                        return index_res;
//                    }
//                    else {
//                        cout << " i need You, my Lord, please help me! not find the key " << key << endl;
//                    }

                        /*
                        auto findf = (bool)buffer->ALEX_DATA_NODE_FLAG_AT(model_res);
                        auto findv = (key_type*)buffer->ALEX_DATA_NODE_PAYLOAD_AT(model_res);
                        auto sortpos = (*buffer)(key) ;
                        if (buffer->slotkey[sortpos] == key) {
                            result_find index_res = result_find(false, buffer->locbitmap[sortpos], buffer,
                                                                sortpos);   // sort_list data
                            return index_res;
                        } else {
                            cout << " i need You, my Lord, please help me! not find the key " << key << endl;
                        }
                        */
                    }
//                    gettimeofday(&ct2, NULL);
//                    ctimeuse = (ct2.tv_sec - ct1.tv_sec) + (double) (ct2.tv_usec - ct1.tv_usec) / 1000000.0;
//                    query_stats->sort_piecetime += ctimeuse;
                    if ( a->delbitmap[resslot]){
                        result_find index_res = result_find(true, a->locbitmap[resslot], a, -1);   // regular data
                        return index_res;
                    }
                    result_find index_res = result_find(true, a->locbitmap[resslot], a, resslot);   // regular data
                    return index_res;
                }

        }

        std::tuple<bool, bool, unsigned long,bool,Leafpiece * > delete_one(const key_type key) {
            // 从 root 开始，即，从root开始的结尾开始

            auto itlevel = innerlevels.end() - 1; // root level 的 iterator

            auto predicted_pos = (*root)(key);
            auto childpieces = (*(--itlevel))->innerpieces;
//                auto search_pos = std::max(exponential_search_upper_bound(childpieces,predicted_pos, key) -1,0);
            auto search_pos = exponential_search_upper_bound(childpieces, predicted_pos, key) - 1;

            // 判断 找到的inner node 是否满足条件
            auto curit = (*(itlevel))->innerpieces.begin() + search_pos;
            Innerpiece *b = *curit;

            for (; --itlevel >= innerlevels.begin();) {
                predicted_pos = (*b)(key);
                childpieces = (*(itlevel))->innerpieces;
                search_pos = exponential_search_upper_bound(childpieces, predicted_pos, key) - 1;

                // 判断 找到的inner node 是否满足条件
                curit = (*(itlevel))->innerpieces.begin() + search_pos;
            }
            // find the leaf level

            b = *curit;
            predicted_pos = (*b)(key);

            search_pos = exponential_search_upper_bound(leaflevel.leafpieces, predicted_pos, key) - 1;

            // 预测位置和最终位置的距离
            auto curleaf = leaflevel.leafpieces.begin() + search_pos;

            Leafpiece *a = *curleaf;
            predicted_pos = (*a)(key);
            // 确定 该leaf 的哪一个链 leaf piece
            while (a->buffer_piece->buffer_piece && key > a->endkey) {
                a = a->buffer_piece;
            }

//            以上的代码定位到m_piece；以下的代码 在m_piece 以及b_piece 中找到要删除的key, 并执行删除

            if (a->buffer_piece == NULL) {    // 这个if 的逻辑是哪种情况？ 直接就在buffer 中find_key: 说明queried_key 大于m_piece 中的所有keys
//                    gettimeofday(&ct1, NULL);
                sortPiece *buffer = (sortPiece *) a;
                auto resslot = predicted_pos * this->buffer_ratio;
                auto model_res = buffer->find_key(key, resslot);
                // 判断要删除的key 是否在该 b_piece 中
                if (model_res != -1){
                    bool in_mem_flag = buffer->ALEX_DATA_NODE_FLAG_AT(model_res);   // 这个是在内存还是在外存的意思； 那判断是有效还是无效key 是用 check_exists
                    auto res = buffer->delete_key(model_res, in_mem_flag);  // 实际返回的是这个值：short_reduced_evict

                    return std::make_tuple(true, in_mem_flag, res,true,buffer);
                }
                else{
                    return std::make_tuple(false, true, 0,true, buffer);
                }
                // 当所有删除的keys 不重复的时候，也就是删除的keys 一定在index 中，则用如下的代码，不需要进行判断
//                bool in_mem_flag = buffer->ALEX_DATA_NODE_FLAG_AT(model_res);   // 这个是在内存还是在外存的意思； 那判断是有效还是无效key 是用 check_exists
//                auto res = buffer->delete_key(key, resslot, in_mem_flag);  // 实际返回的是这个值：short_reduced_evict
//                return std::make_tuple(in_mem_flag, res);
            }
            else {
                predicted_pos = (*a)(key);

                int slotlo = PGM_SUB_EPS(predicted_pos, ErrorRecursive + 1);
                slotlo = std::min<int>(slotlo,(a->slotkey.size() - 1));
                int slothi = PGM_ADD_EPS(predicted_pos, ErrorRecursive, (a->slotkey.size() - 1));
                int resslot = slotlo;
                for (; resslot < slothi; resslot++) {
                    if (key <= a->slotkey[resslot]) {
                        break;
                    }
                }
//                    gettimeofday(&ct1, NULL);
                if (key != a->slotkey[resslot]) {
                    // find in sort_list
                    while (a->buffer_piece->buffer_piece != NULL) {
                        a = a->buffer_piece;
                    }
                    predicted_pos = (*(Leafpiece*)(*curleaf))(key);
                    sortPiece *buffer = (sortPiece *) a->buffer_piece;
                    resslot = predicted_pos * this->buffer_ratio;
                    auto model_res = buffer->find_key(key, resslot);

                    if (model_res != -1){
                        bool in_mem_flag = buffer->ALEX_DATA_NODE_FLAG_AT(model_res);   // 这个是在内存还是在外存的意思； 那判断是有效还是无效key 是用 check_exists
                        auto res = buffer->delete_key(model_res, in_mem_flag);  // 实际返回的是这个值：short_reduced_evict
                        return std::make_tuple(true, in_mem_flag, res,true, buffer);
                    }
                    else{
                        return std::make_tuple(false, true, 0,true, buffer);
                    }

//                    bool in_mem_flag = buffer->ALEX_DATA_NODE_FLAG_AT(model_res);
//                    auto res_pos = buffer->delete_key(key, resslot, in_mem_flag);  // short_reduced_evict
//                    return std::make_tuple(in_mem_flag, res_pos);

                }
                // 说明要删除的key 是在m_piece 中
                if (a->delbitmap[resslot]){  // 为true, 说明已经删除了
                    return std::make_tuple(false, true, 0, true, a);
                }
                else{
                    bool in_mem_flag = a->locbitmap[resslot];
                    auto [m_size, res] = a->delete_key(resslot,in_mem_flag);   // short_reduced_evict
                    return std::make_tuple(true, in_mem_flag, res, m_size, a);

                }
//                bool in_mem_flag = a->locbitmap[resslot];
//                auto res = a->delete_key(key,resslot,in_mem_flag);   // short_reduced_evict
//                return std::make_tuple(in_mem_flag, res);
            }

        }

        forceinline range_key_result_find search_one_lower(key_type key) {
            // 找到 第一个大于或者等于 search_key 的piece (包括leaf，buffer)， 及对应的位置
//            struct timeval ct1, ct2;
//            double ctimeuse;
//            if (vstats.innernum > 1) {
                auto itlevel = innerlevels.end() - 1; // root level 的 iterator
                auto predicted_pos = (*root)(key);
                auto childpieces = (*(--itlevel))->innerpieces;
                auto search_pos = exponential_search_upper_bound(childpieces, predicted_pos, key) - 1;
//                if (search_pos <0)
//                    search_pos = 0;
                // 判断 找到的inner node 是否满足条件
                auto curit = (*(itlevel))->innerpieces.begin() + search_pos;
                Innerpiece *b = *curit;
//                if (b->startkey > key )
//                    cout << "i need You, my lovely Lord, b->startkey > key " << endl;

                for (; --itlevel >= innerlevels.begin();) {
                    predicted_pos = (*b)(key);
                    childpieces = (*(itlevel))->innerpieces;
                    search_pos = exponential_search_upper_bound(childpieces, predicted_pos, key) - 1;
//                    if (search_pos <0)
//                        search_pos = 0;
                    curit = (*(itlevel))->innerpieces.begin() + search_pos;
                }

                b = *curit;
                predicted_pos = (*b)(key);
                search_pos = exponential_search_upper_bound(leaflevel.leafpieces, predicted_pos, key) - 1;
//                if (search_pos <0)
//                    search_pos = 0;
                auto curleaf = leaflevel.leafpieces.begin() + search_pos;

                Leafpiece *a = *curleaf;
                predicted_pos = (*a)(key);
                // 确定 该leaf 的哪一个链 leaf piece
                while (key > a->endkey && a->buffer_piece) {
                    a = a->buffer_piece;
                }
                if (a->buffer_piece == NULL) {
                    // 这种情况，说明，leaf piece 中所有的key 都小于search key
//                    将 Leafpiece* leaf_a = *curleaf 给findleafs，并将slot = piece.slotkey.size 给 slots;
                    Leafpiece *leaf_a = *curleaf;
                    range_key_result_find index_res = range_key_result_find(search_pos);
                    index_res.findleafs.emplace_back(leaf_a);
                    index_res.slots.push_back(leaf_a->slotkey.size());
                    index_res.findleafs.emplace_back(a);

                    sortPiece *buffer = (sortPiece *) a;
                    auto resslot = predicted_pos * this->buffer_ratio;

                    auto model_res = buffer->find_key_lower(key, resslot);
//                        auto model_res = buffer->find_key(key, resslot);  // 精确找到key，不适合用于 range query
//                        if (model_res == buffer->data_capacity_)
//                            model_res -= 1;
                    index_res.slots.push_back(model_res);
                    return index_res;
                }
                else {
                    predicted_pos = (*a)(key);

                    int slotlo = PGM_SUB_EPS(predicted_pos, ErrorRecursive + 1);
                    slotlo = std::min<int>(slotlo,(a->slotkey.size() - 1));
                    int slothi = PGM_ADD_EPS(predicted_pos, ErrorRecursive, (a->slotkey.size() - 1));
                    int resslot = slotlo;
                    for (; resslot < slothi; resslot++) {
                        if (key <= a->slotkey[resslot]) {
                            break;
                        }
                    }
                    range_key_result_find index_res = range_key_result_find(search_pos);
                    index_res.findleafs.emplace_back(a);
                    index_res.slots.push_back(resslot);

                    // find in sort_list
                    while (a->buffer_piece->buffer_piece != NULL) {
                        a = a->buffer_piece;
                    }

                    index_res.findleafs.emplace_back(a->buffer_piece);
                    sortPiece *buffer = (sortPiece *) (sortPiece *) a->buffer_piece;
                    resslot = predicted_pos * this->buffer_ratio;


                    if (buffer->num_keys_ ==0){
                        auto model_res = -1;
                        index_res.slots.push_back(model_res);
                        return index_res;
                    }

                    else{
                        auto model_res = buffer->find_key_lower(key, resslot)+1;  // find_key_lower 似乎是找到第一个小于等于key 的位置
//                        auto model_res = buffer->find_key(key, resslot);  // 精确找到key，不适合用于 range query
//                        if (model_res == buffer->data_capacity_)
//                            model_res -= 1;
                        index_res.slots.push_back(model_res);
                        return index_res;
                    }


                    }

            }

        forceinline range_key_result_find search_one_scan_lower(key_type key) {
            // 找到 第一个大于或者等于 search_key 的piece (包括leaf，buffer)， 及对应的位置
//            struct timeval ct1, ct2;
//            double ctimeuse;
//            if (vstats.innernum > 1) {
                auto itlevel = innerlevels.end() - 1; // root level 的 iterator
                auto predicted_pos = (*root)(key);
                auto childpieces = (*(--itlevel))->innerpieces;
                auto search_pos = exponential_search_upper_bound(childpieces, predicted_pos, key) - 1;
//                if (search_pos <0)
//                    search_pos = 0;
                // 判断 找到的inner node 是否满足条件
                auto curit = (*(itlevel))->innerpieces.begin() + search_pos;
                Innerpiece *b = *curit;
//                if (b->startkey > key )
//                    cout << "i need You, my lovely Lord, b->startkey > key " << endl;

                for (; --itlevel >= innerlevels.begin();) {
                    predicted_pos = (*b)(key);
                    childpieces = (*(itlevel))->innerpieces;
                    search_pos = exponential_search_upper_bound(childpieces, predicted_pos, key) - 1;
//                    if (search_pos <0)
//                        search_pos = 0;
                    curit = (*(itlevel))->innerpieces.begin() + search_pos;
                }

                b = *curit;
                predicted_pos = (*b)(key);
                search_pos = exponential_search_upper_bound(leaflevel.leafpieces, predicted_pos, key) - 1;
//                if (search_pos <0)
//                    search_pos = 0;
                auto curleaf = leaflevel.leafpieces.begin() + search_pos;

                Leafpiece *a = *curleaf;
                predicted_pos = (*a)(key);
                // 确定 该leaf 的哪一个链 leaf piece
                while (key > a->endkey && a->buffer_piece) {
                    a = a->buffer_piece;
                }
                if (a->buffer_piece == NULL) {
                    // 这种情况，说明，leaf piece 中所有的key 都小于search key
//                    将 Leafpiece* leaf_a = *curleaf 给findleafs，并将slot = piece.slotkey.size 给 slots;
                    Leafpiece *leaf_a = *curleaf;
                    range_key_result_find index_res = range_key_result_find(search_pos);
                    index_res.findleafs.emplace_back(leaf_a);
                    index_res.slots.push_back(leaf_a->slotkey.size());
                    index_res.findleafs.emplace_back(a);

                    sortPiece *buffer = (sortPiece *) a;
                    auto resslot = predicted_pos * this->buffer_ratio;

//                    auto model_res = buffer->find_key_lower(key, resslot);
                    auto model_res = buffer->find_key_upper(key, resslot);

                    index_res.slots.push_back(model_res);
                    return index_res;
                }
                else {
                    predicted_pos = (*a)(key);

                    int slotlo = PGM_SUB_EPS(predicted_pos, ErrorRecursive + 1);
                    slotlo = std::min<int>(slotlo,(a->slotkey.size() - 1));
                    int slothi = PGM_ADD_EPS(predicted_pos, ErrorRecursive, (a->slotkey.size() - 1));
                    int resslot = slotlo;
                    for (; resslot < slothi; resslot++) {
                        if (key <= a->slotkey[resslot]) {
                            break;
                        }
                    }
                    range_key_result_find index_res = range_key_result_find(search_pos);
                    index_res.findleafs.emplace_back(a);
                    index_res.slots.push_back(resslot);

                    // find in sort_list
                    while (a->buffer_piece->buffer_piece != NULL) {
                        a = a->buffer_piece;
                    }

                    index_res.findleafs.emplace_back(a->buffer_piece);
                    sortPiece *buffer = (sortPiece *) (sortPiece *) a->buffer_piece;
                    resslot = predicted_pos * this->buffer_ratio;

                    if (buffer->num_keys_ ==0){
//                        auto model_res = -1;
                        index_res.slots.push_back(buffer->data_capacity_);
                        return index_res;
                    }
                    else{
//                        auto model_res = buffer->find_key_lower(key, resslot);
                        auto model_res = buffer->find_key_upper(key, resslot);  // 找到大于等于key 的第一个有效位置
//                        if (model_res == buffer->data_capacity_)
//                            model_res -= 1;
                        if (model_res == -1)
                            model_res = buffer->find_first_1();
                        index_res.slots.push_back(model_res);
//                        auto mkey = buffer->ALEX_DATA_NODE_KEY_AT(model_res);
//                        if (mkey < key)
//                            COUT_THIS("thank You, my lovely Lord");
                        return index_res;
                    }

                }
//            }
//            else{
//
//                // 只有一层leaf
////                static constexpr size_t linear_search_threshold = 8 * 64 / sizeof(Innerpiece);   //sizeof(Innerpiece)
//
//                auto predicted_pos = (*root)(key);
//
//                auto search_pos = exponential_search_upper_bound(leaflevel.leafpieces, predicted_pos, key) - 1;
////                if (search_pos <0)
////                    search_pos = 0;
//                auto curleaf = leaflevel.leafpieces.begin() + search_pos;   // leaf level
//
//                Leafpiece *a = *curleaf;
//
//                // 确定 该leaf 的哪一个链 leaf piece
//                while (key > a->endkey && a->buffer_piece) {
//                    a = a->buffer_piece;
//                }
//                if (a->buffer_piece == NULL) {
//                    // 这种情况，说明，leaf piece 中所有的key 都小于search key
////                    将 Leafpiece* leaf_a = *curleaf 给findleafs，并将slot = piece.slotkey.size 给 slots;
//                    Leafpiece *leaf_a = *curleaf;
//                    range_key_result_find index_res = range_key_result_find(search_pos);
//                    index_res.findleafs.emplace_back(leaf_a);
//                    index_res.slots.push_back(leaf_a->slotkey.size() - 1);
//                    index_res.findleafs.emplace_back(a);
//                    sortPiece *buffer = (sortPiece *) a;
//                    auto resslot = predicted_pos * this->buffer_ratio;
//
//                    auto model_res = buffer->find_key(key, resslot);
//                    index_res.slots.push_back(model_res);
//                    return index_res;
//                }
//                else {
//                    predicted_pos = (*a)(key);
//
//                    int slotlo = PGM_SUB_EPS(predicted_pos, ErrorRecursive + 1);
//                    slotlo = std::min<int>(slotlo,(a->slotkey.size() - 1));
//                    int slothi = PGM_ADD_EPS(predicted_pos, ErrorRecursive, (a->slotkey.size() - 1));
//                    int resslot = slotlo;
//                    for (; resslot < slothi; resslot++) {
//                        if (key <= a->slotkey[resslot]) {
//                            break;
//                        }
//                    }
//                    range_key_result_find index_res = range_key_result_find(search_pos);
//                    index_res.findleafs.emplace_back(a);
//                    index_res.slots.push_back(resslot);
//
//                    // find in sort_list
//                    while (a->buffer_piece->buffer_piece != NULL) {
//                        a = a->buffer_piece;
//                    }
//
//                    index_res.findleafs.emplace_back(a->buffer_piece);
//                    sortPiece *buffer = (sortPiece *) (sortPiece *) a->buffer_piece;
//                    resslot = predicted_pos * this->buffer_ratio;
//
//                    if (buffer->num_keys_ ==0){
////                        auto model_res = -1;
//                        index_res.slots.push_back(buffer->data_capacity_);
//                        return index_res;
//                    }
//                    else{
//                        auto model_res = buffer->find_key_lower(key, resslot);
////                        if (model_res == buffer->data_capacity_)
////                            model_res -= 1;
//                        index_res.slots.push_back(model_res);
//                        return index_res;
//                    }
//                }
//
//            }
        }

        /* * for large error, that may have only one leaf piece
         * /
//        forceinline result_find search_one_leaf(const key_type key) {
//            // 从 root 开始，即，从root开始的结尾开始
//            struct timeval ct1, ct2;
//            double ctimeuse;
//            if (vstats.innernum > 1){
//                auto itlevel = innerlevels.end() -1; // root level 的 iterator
//                // internal level use exponential search, the parameters used in exponential_search are: predict_pos, search key
//                /*
//                 * 1. get the predicted pos
//                 * 2. get the [l,r] used in binary search or, get the [l,r]
//                 * 3. use the binary search in [l,r] to get the final position
//                 */
//
//                auto predicted_pos = (*root)(key);
//                auto childpieces = (*(--itlevel))->innerpieces;
//                auto search_pos = exponential_search_upper_bound(childpieces,predicted_pos, key) -1;
//                if (search_pos <0)
//                    search_pos = 0;
//                // 判断 找到的inner node 是否满足条件
//                auto curit = (*(itlevel))->innerpieces.begin() + search_pos;
//                Innerpiece* b = *curit;
////                if (b->startkey > key )
////                    cout << "i need You, my lovely Lord, b->startkey > key " << endl;
//
//                for (; --itlevel >= innerlevels.begin();){
//                    predicted_pos = (*b)(key);
//                    childpieces = (*(itlevel))->innerpieces;
//                    search_pos = exponential_search_upper_bound(childpieces,predicted_pos, key) -1;
//                    if (search_pos <0)
//                        search_pos = 0;
//                    // 判断 找到的inner node 是否满足条件
//                    curit = (*(itlevel))->innerpieces.begin() + search_pos;
//                }
//
//                // find the leaf level
////                auto curleaf = leaflevel.leafpieces.begin();
//
//                b = *curit;
//                predicted_pos = (*b)(key);
//
//                search_pos = exponential_search_upper_bound(leaflevel.leafpieces,predicted_pos, key) -1;
//                if (search_pos <0)
//                    search_pos = 0;
//                auto curleaf = leaflevel.leafpieces.begin() + search_pos;
//
//                Leafpiece* a = *curleaf;
//
//                result_find index_res= result_find(false,a->locbitmap[0],a,0);   // not find, then, insert into bufferpiece,
//                index_res.leafslot = search_pos;
//
//                return index_res;
//            }
//            else{
//
//                // 只有一层leaf
//                static constexpr size_t linear_search_threshold = 8 * 64 / sizeof(Innerpiece);   //sizeof(Innerpiece
//
//                auto predicted_pos = (*root)(key);
//
//
//                auto search_pos = exponential_search_upper_bound(leaflevel.leafpieces,predicted_pos, key) -1;
//                if (search_pos <0)
//                    search_pos = 0;
//                auto curleaf = leaflevel.leafpieces.begin() + search_pos;  // leaf level
//
//                Leafpiece* a = *curleaf;
//
//                result_find index_res= result_find(false,a->locbitmap[0],a,0);   // 默认 数据是不重复的，
//                // 如果数据是允许重复的，这里需要加判断条件，如果找到了，进行update，if not find, then, insert into bufferpiece,
//                index_res.leafslot = search_pos;
//
//                return index_res;
//            }
//        }


        forceinline std::pair<Leafpiece*,int> search_one_leaf(const key_type key) {
            // 从 root 开始，即，从root开始的结尾开始
//            if (vstats.innernum > 1) {
                auto itlevel = innerlevels.end() -1; // root level 的 iterator
                // internal level use exponential search, the parameters used in exponential_search are: predict_pos, search key
                /*
                 * 1. get the predicted pos
                 * 2. get the [l,r] used in binary search or, get the [l,r]
                 * 3. use the binary search in [l,r] to get the final position
                 */

                auto predicted_pos = (*root)(key);
                auto childpieces = (*(--itlevel))->innerpieces;
//                auto search_pos = std::max(exponential_search_upper_bound(childpieces,predicted_pos, key) -1,0) ;
            auto search_pos = exponential_search_upper_bound(childpieces,predicted_pos, key) -1 ;
//            if (search_pos <0)
//                search_pos = 0;
                // 判断 找到的inner node 是否满足条件
                auto curit = (*(itlevel))->innerpieces.begin() + search_pos;
                Innerpiece* b = *curit;
//                if (b->startkey > key )
//                    cout << "i need You, my lovely Lord, b->startkey > key " << endl;

                for (; --itlevel >= innerlevels.begin();){
                    predicted_pos = (*b)(key);
                    childpieces = (*(itlevel))->innerpieces;
                    search_pos = std::max(exponential_search_upper_bound(childpieces,predicted_pos, key) -1,0);
//                if (search_pos <0)
//                    search_pos = 0;
                    // 判断 找到的inner node 是否满足条件
                    curit = (*(itlevel))->innerpieces.begin() + search_pos;
                }

                // find the leaf level
//                auto curleaf = leaflevel.leafpieces.begin();

                b = *curit;
                predicted_pos = (*b)(key);

                search_pos = std::max(exponential_search_upper_bound(leaflevel.leafpieces,predicted_pos, key) -1,0);
//            if (search_pos <0)
//                search_pos = 0;
//            auto curleaf = leaflevel.leafpieces.begin() + search_pos;

                Leafpiece* a = *(leaflevel.leafpieces.begin() + search_pos);

//            result_find index_res= result_find(false,a->locbitmap[0],a,0);   // not find, then, insert into bufferpiece,
//            index_res.leafslot = search_pos;

                return {a,search_pos};
//            }
//            else{
//                auto predicted_pos = (*root)(key);
//
//                auto search_pos = exponential_search_upper_bound(leaflevel.leafpieces, predicted_pos, key) - 1;
////                if (search_pos <0)
////                    search_pos = 0;
//                Leafpiece* a = *(leaflevel.leafpieces.begin() + search_pos);
//
////            result_find index_res= result_find(false,a->locbitmap[0],a,0);   // not find, then, insert into bufferpiece,
////            index_res.leafslot = search_pos;
//
//                return {a,search_pos};
//            }

//            auto itlevel = innerlevels.end() -1; // root level 的 iterator
//            // internal level use exponential search, the parameters used in exponential_search are: predict_pos, search key
//            /*
//             * 1. get the predicted pos
//             * 2. get the [l,r] used in binary search or, get the [l,r]
//             * 3. use the binary search in [l,r] to get the final position
//             */
//
//            auto predicted_pos = (*root)(key);
//            auto childpieces = (*(--itlevel))->innerpieces;
//            auto search_pos = std::max(exponential_search_upper_bound(childpieces,predicted_pos, key) -1,0) ;
////            if (search_pos <0)
////                search_pos = 0;
//            // 判断 找到的inner node 是否满足条件
//            auto curit = (*(itlevel))->innerpieces.begin() + search_pos;
//            Innerpiece* b = *curit;
////                if (b->startkey > key )
////                    cout << "i need You, my lovely Lord, b->startkey > key " << endl;
//
//            for (; --itlevel >= innerlevels.begin();){
//                predicted_pos = (*b)(key);
//                childpieces = (*(itlevel))->innerpieces;
//                search_pos = std::max(exponential_search_upper_bound(childpieces,predicted_pos, key) -1,0);
////                if (search_pos <0)
////                    search_pos = 0;
//                // 判断 找到的inner node 是否满足条件
//                curit = (*(itlevel))->innerpieces.begin() + search_pos;
//            }
//
//            // find the leaf level
////                auto curleaf = leaflevel.leafpieces.begin();
//
//            b = *curit;
//            predicted_pos = (*b)(key);
//
//            search_pos = std::max(exponential_search_upper_bound(leaflevel.leafpieces,predicted_pos, key) -1,0);
////            if (search_pos <0)
////                search_pos = 0;
////            auto curleaf = leaflevel.leafpieces.begin() + search_pos;
//
//            Leafpiece* a = *(leaflevel.leafpieces.begin() + search_pos);
//
////            result_find index_res= result_find(false,a->locbitmap[0],a,0);   // not find, then, insert into bufferpiece,
////            index_res.leafslot = search_pos;
//
//            return {a,search_pos};
        }


        htaprange_res_find search_range(key_type firstkey, key_type lastkey) {
            struct timeval tqt1, tqt2, xt1, xt2, plt1, plt2, prdt1, prdt2;
            double pltimeuse, prdtimeuse, pwdtimeuse, timeuse, xtimeuse;  // 读磁盘的时间，写磁盘的时间
            pair<key_type, key_type *> res;
            gettimeofday(&xt1, NULL);
            auto index_res = search_one_lower(firstkey);
            gettimeofday(&xt2, NULL);
            // according to the index_res of first key to define the laskkey's index_res

            htaprange_res_find rangeresult = htaprange_res_find(firstkey, lastkey, index_res.slots[0],index_res.slots[1]);   // what's the meaning of the first parameter?
            // 根据first key的位置定位 lastkey 的位置
            // 找到lastkey 所在leaf的位置，2 找到lastkey 所在buffer 的位置。 都是返回最后一个小于或等于 lastkey 的position
            auto curleaf = leaflevel.leafpieces.begin()+index_res.leafslot;
            Leafpiece *a = index_res.findleafs[0];
            auto leaf_res = a->find_upper_in_leaf(lastkey, Error);
            while (leaf_res.first == false && a!=this->m_tailleaf)  {

                if (a->buffer_piece->buffer_piece!= NULL){
                    rangeresult.findleafs.push_back((a));
                    rangeresult.leafslots.push_back(index_res.leafslot);
                    a = a->buffer_piece;
                }
                else{
                    rangeresult.findleafs.push_back((a));
                    rangeresult.leafslots.push_back(index_res.leafslot);
                    index_res.leafslot ++;
                    a = *(leaflevel.leafpieces.begin()+index_res.leafslot);
                }

                leaf_res = a->find_upper_in_leaf(lastkey, Error);
            }
            if (a->buffer_piece!=NULL)
                rangeresult.findleafs.push_back(a);
            // 这里不加区分，就放入curleaf， 不合理，因为findsleafs[0] 有可能不是leaflevel slot 上的第一个leaf piece
            rangeresult.leafslots.push_back(index_res.leafslot++);
            rangeresult.slotleaf1 = leaf_res.second;
            int leaf_num= rangeresult.findleafs.size();
            if (leaf_num>1){
                sortPiece * buffer = (sortPiece *) (index_res.findleafs[1]);
                if (buffer->num_keys_ >0)
                    rangeresult.findbuffers.emplace_back(buffer);

                // 检查新的 leaf 的buffer
                for (int i = 1; i<leaf_num;i++){
                    buffer = (sortPiece *) (rangeresult.findleafs[i]->buffer_piece);
                    while (buffer->buffer_piece != NULL){
                        Leafpiece * curleaf = (Leafpiece *)(buffer);
                        buffer = (sortPiece *) (curleaf->buffer_piece);
                    }
                    if (buffer->num_keys_ > 0){   // = 0 说明buffer 中没有数据   446181678460  446345317732  446346181696
                        auto predict_pos = (*rangeresult.findleafs[i])(lastkey)*buffer_ratio;
                        auto buffer_res = buffer->find_upper_in_leaf(lastkey, predict_pos);
                        rangeresult.findbuffers.emplace_back(buffer);
                        rangeresult.slotbuffer1 = buffer_res;
                    }
                    else if(rangeresult.findbuffers.size()){
                        rangeresult.slotbuffer1 = rangeresult.findbuffers.back()->data_capacity_;
                    }

                }
            }
            else{
                sortPiece *buffer = (sortPiece *) (index_res.findleafs[1]);
                if (buffer->num_keys_ > 0){   // = 0 说明buffer 中没有数据
//                    auto mkey1 = buffer->ALEX_DATA_NODE_KEY_AT(331);
//                    auto mkey2 = buffer->ALEX_DATA_NODE_KEY_AT(332);
                    auto predict_pos =  (*rangeresult.findleafs[0])(lastkey)*buffer_ratio;// 计算 predict_pos;
                    auto buffer_res = buffer->find_upper_in_leaf(lastkey, predict_pos);
                    rangeresult.findbuffers.emplace_back(buffer);
                    rangeresult.slotbuffer1 = buffer_res;
                }
            }


//            cout << "i trust in You, my Lord" << endl;
            return rangeresult;

        }
//
        htaprange_res_find search_scan(key_type firstkey, const int scan_num) {
            struct timeval tqt1, tqt2, xt1, xt2, plt1, plt2, prdt1, prdt2;
            double pltimeuse, prdtimeuse, pwdtimeuse, timeuse, xtimeuse;  // 读磁盘的时间，写磁盘的时间
            pair<key_type, key_type *> res;
//            gettimeofday(&xt1, NULL);
            auto index_res = search_one_scan_lower(firstkey);
//            gettimeofday(&xt2, NULL);
            // according to the index_res of first key to define the laskkey's index_res
            int remaining = scan_num;
            htaprange_res_find rangeresult = htaprange_res_find(firstkey, index_res.slots[0], index_res.slots[1]);   // what's the meaning of the first parameter?
            // 根据first key的位置定位 lastkey 的位置
            // 找到lastkey 所在leaf的位置，2 找到lastkey 所在buffer 的位置。 都是返回最后一个小于或等于 lastkey 的position
//            auto curleaf = leaflevel.leafpieces.begin()+index_res.leafslot;

            // do scan
            // model piece
//            Leafpiece *a = index_res.findleafs[0];
            Leafpiece *findpiece = index_res.findleafs[0];
            auto start_leaf = index_res.slots[0];
            int leafmax = findpiece->slotkey.size();

            // buffer
            sortPiece * bufferpiece = (sortPiece *) (index_res.findleafs[1]);
//            int bufmax = bufferpiece->slotkey.size()-1;
            auto start_buffer = index_res.slots[1];
//            auto start_buffer1 = index_res.slots[1];
//            if (start_buffer <0 ){
//                start_buffer1 ++;
//                start_buffer++;
//            }

            // 判断是否移动到下一个leaf
            if (start_leaf == findpiece->slotkey.size()  ){
//                rangeresult.findleafs.push_back(findpiece);
                rangeresult.slotleaf0 = 0;
                // 下一个leafslot 还是 下一个model piece
                if (findpiece->buffer_piece == bufferpiece){
                    // next leaf slot
                    findpiece = *(leaflevel.leafpieces.begin()+ (++index_res.leafslot));
//                    while( start_buffer1< bufferpiece->data_capacity_  && remaining ){
//                        if (bufferpiece->check_exists(start_buffer1++) ){
////                            auto mkey = bufferpiece->ALEX_DATA_NODE_KEY_AT(start_buffer1-1);
//                            COUT_THIS("thank You, my Lord, please come, Holy Spirit!");
//                        }
//                    }
                    if (start_buffer < bufferpiece->data_capacity_){
                        while( start_buffer< bufferpiece->data_capacity_  && remaining ){
                            if (bufferpiece->check_exists(start_buffer) && bufferpiece->ALEX_DATA_NODE_KEY_AT(start_buffer)>=firstkey ){
                                remaining --;
                            }
                            start_buffer ++;
                        }
                        rangeresult.findbuffers.emplace_back(bufferpiece);
                        bufferpiece = (sortPiece *) (findpiece->buffer_piece);
                        while (bufferpiece->buffer_piece != NULL)
                            bufferpiece = (sortPiece *) bufferpiece->buffer_piece;
                        start_leaf = 0;
                        // 找到该buffer 中的第一个有效record 的位置
                        start_buffer = bufferpiece->find_first_1();
//                        rangeresult.slotbuffer0 = start_buffer;
//                        start_buffer1 = 0;
                    }
                    else{
                        bufferpiece = (sortPiece *) (findpiece->buffer_piece);
                        while (bufferpiece->buffer_piece != NULL)
                            bufferpiece = (sortPiece *) bufferpiece->buffer_piece;
                        start_leaf = 0;
                        // 找到该buffer 中的第一个有效record 的位置
                        start_buffer = bufferpiece->find_first_1();
                        rangeresult.slotbuffer0 = start_buffer;
//                        start_buffer1 = 0;
                    }

//                        if (start_buffer != start_buffer1){
//                            COUT_THIS("thank You, my lovely Lord");
//                        }
                }
                else{ // next model piece
                    rangeresult.findleafs.push_back((findpiece));
                    rangeresult.leafslots.push_back(index_res.leafslot);
                    findpiece = findpiece->buffer_piece;
                    start_leaf = 0;
                }
            }


            while(remaining>0 && start_leaf< findpiece->slotkey.size() || remaining>0 && start_buffer < bufferpiece->data_capacity_){


                if (findpiece->slotkey[start_leaf] < bufferpiece->key_slots_[start_buffer] && start_leaf< findpiece->slotkey.size()
                || start_buffer > (bufferpiece->data_capacity_-1)){
                    start_leaf ++;
                    remaining --;
                }

                else{
//                    if (start_buffer != start_buffer1)
//                        COUT_THIS("my Lord, please be with Brother Mingde");
                    while(start_buffer < bufferpiece->data_capacity_  && !bufferpiece->check_exists(start_buffer++) ){
                    }
//                    while(start_buffer1< bufferpiece->data_capacity_  && !bufferpiece->check_exists(start_buffer1++) ){
//                    }
//                    if (start_buffer != start_buffer1)
//                        COUT_THIS("my Lord, please be with Brother Mingde");
//                    start_buffer ++;
                    remaining --;
                }


                // 判断是否移动到下一个leaf
                if (start_leaf == findpiece->slotkey.size()  ){

//                    rangeresult.findbuffers.push_back(bufferpiece);
                    // 下一个leafslot 还是 下一个model piece
                    if (findpiece->buffer_piece == bufferpiece){
                        // next leaf slot
                        if (++index_res.leafslot==vstats.leaves){
                            remaining = 0;
                            start_leaf --;
                            break;
                        }
                        rangeresult.findleafs.push_back(findpiece);
                        findpiece = *(leaflevel.leafpieces.begin()+ index_res.leafslot);
//                        while( start_buffer1< bufferpiece->data_capacity_  && remaining ){
//                            if (bufferpiece->check_exists(start_buffer1++) )
//                                auto m = remaining;
//                        }
                        while( start_buffer< bufferpiece->data_capacity_  && remaining ){
                            if (bufferpiece->check_exists(start_buffer++) )
                                remaining --;
                        }
                        rangeresult.findbuffers.emplace_back(bufferpiece);
                        bufferpiece = (sortPiece *) (findpiece->buffer_piece);
                        while (bufferpiece->buffer_piece != NULL)
                            bufferpiece = (sortPiece *) bufferpiece->buffer_piece;
                        start_leaf = 0;
                        // 找到该buffer 中的第一个有效record 的位置
                        start_buffer = bufferpiece->find_first_1();
//                        start_buffer1 = 0;
//                        if (start_buffer != start_buffer1){
//                            COUT_THIS("thank You, my lovely Lord");
//                        }
                    }
                    else{ // next model piece
                        rangeresult.findleafs.push_back((findpiece));
                        rangeresult.leafslots.push_back(index_res.leafslot);
                        findpiece = findpiece->buffer_piece;
                        start_leaf = 0;
                    }
                }
            }

            rangeresult.findleafs.push_back(findpiece);
            rangeresult.findbuffers.push_back(bufferpiece);
            rangeresult.slotleaf1 = start_leaf;
            rangeresult.slotbuffer1 = start_buffer;
            rangeresult.leafslots.push_back(index_res.leafslot);

//            cout << "i trust in You, my Lord" << endl;
            return rangeresult;

        }

//
//        template<class lru_type>
//        forceinline bool insert_one(const key_type key, key_type* payload,lru_type interchain)  {
//            // 判断一下，如果 小于 last leaf 的end key, 则插入到sortpiece 中
//            if (key < this->m_tailleaf->endkey){
//                // 找到属于哪个leaf piece
//                // 插入到 该leaf piece 的buffer piece 中，
//                // 判断插入到leaf piece 中，是否超出阈值
//                auto search_res = this->search_one_leaf(key);
//                Leafpiece* a = search_res.findleaf;
//                int link_num = 0;
//                while (a->buffer_piece->buffer_piece != NULL){
//                    a = a->buffer_piece;
//                    link_num ++;
//                }
////                if ( a->startkey == 11599762525505 ){   // 16665361923015
////                    COUT_THIS("thank You, my lovely Lord!");
////                    }
//
//                sortPiece* bufferpiece = (sortPiece*) (a->buffer_piece);
//
//                auto lrukey = bufferpiece->key_slots_[0];
//                auto predict_pos = (*search_res.findleaf)(key) * this->buffer_ratio;  // normalize to buffer size;
//                //std::pair<int, int> positions = bufferpiece->find_insert_position(key, predict_pos);
//                auto mflags = bufferpiece->model_insert(key,payload,predict_pos);
//                if (mflags.first)
//                {  // findleaf rebuild
////                    if ( a->startkey == 14790530759750 ){   // 16665361923015
////                        COUT_THIS("thank You, my lovely Lord!");
////                    }
//                    if (mflags.second){
//                        interchain->remove( lrukey);
////                        cout << "i need You, my Lord" <<endl;
//                    }
//                    search_res.findleaf->mlearnedMergeSort(ErrorRecursive, interchain);
//                    // retrain model, get the new slope and intercept
//                    int newmodels = internal::make_segmentation( this,search_res.findleaf,Error,search_res.findleaf->slotkey,interchain,search_res.leafslot) - link_num;
//                    this->newmodels += newmodels;
//                    this->retrain ++;
//                }
//                else{
//                    if (mflags.second){
//                        interchain->remove( lrukey);
//                    }
//                    auto newlrukey = bufferpiece->key_slots_[0];
////                    if (key != newlrukey)
////                        COUT_THIS("please come, my lovely Lord, key != newlrukey")
//                    interchain->put(newlrukey,bufferpiece);
//                }
////                if ( a->startkey == 13277990357497 ){   // 16665361923015
////
////                    COUT_THIS("Jesus, You are with me");
////                    auto key1 = bufferpiece->ALEX_DATA_NODE_KEY_AT(1);
////                    auto key2 = bufferpiece->ALEX_DATA_NODE_KEY_AT(2);
////                    auto key3 = bufferpiece->ALEX_DATA_NODE_KEY_AT(59);
////                    auto key4 = bufferpiece->ALEX_DATA_NODE_KEY_AT(802);
////                    auto key5 = bufferpiece->ALEX_DATA_NODE_KEY_AT(803);
////                    auto key6 = bufferpiece->ALEX_DATA_NODE_KEY_AT(804);
////                    auto key7 = bufferpiece->ALEX_DATA_NODE_KEY_AT(805);
////                    auto key8 = bufferpiece->ALEX_DATA_NODE_KEY_AT(806);
////                    auto key9 = bufferpiece->ALEX_DATA_NODE_KEY_AT(807);
////                    auto key10 = bufferpiece->ALEX_DATA_NODE_KEY_AT(808);
////                    auto key11 = bufferpiece->ALEX_DATA_NODE_KEY_AT(507);
////                    auto key12 = bufferpiece->ALEX_DATA_NODE_KEY_AT(558);
////                    cout << "O Lord, teach me where and how to seek You" <<endl;
////                }
//            }
//            else{
//                this->append_one(key, payload,Error, interchain);
//            }
//        }


        template<class lru_type>
        forceinline bool insert(const key_type key, key_type* payload,lru_type interchain)  {

            // 找到属于哪个leaf piece
            // 插入到 该leaf piece 的buffer piece 中，
            // 判断插入到leaf piece 中，是否超出阈值
            auto search_res = this->search_one_leaf(key);
            Leafpiece* a = search_res.first;
            int link_num = 0;
//            while (a->buffer_piece->buffer_piece != NULL){
//                a = a->buffer_piece;
//                link_num ++;
//            }

            sortPiece* bufferpiece = (sortPiece*) (a->buffer_piece);
            while(bufferpiece->buffer_piece!=NULL){
                bufferpiece =  (sortPiece*)bufferpiece->buffer_piece;
//                link_num ++;
            }

            // 判断是否进行 single piece retraining
            if (bufferpiece -> num_keys_ == bufferpiece->data_capacity_){
//                auto lrukey = bufferpiece->key_slots_[0];
                interchain->remove( bufferpiece->key_slots_[0]);
                std::vector<key_type> retrain_keys;
                std::vector<void*> retrain_datas;
                std::vector<bool> retrain_flags;
                a->mlearnedMergeSort(ErrorRecursive, interchain,retrain_keys,retrain_datas,retrain_flags);
                // retrain model, get the new slope and intercept
                bufferpiece = internal::make_segmentation( this,Error,retrain_keys,retrain_flags,retrain_datas,interchain,search_res.second) ;
//                this->newmodels += (size_buffer.first- link_num);
                this->retrain ++;
                // 执行插入
                a = leaflevel.leafpieces[search_res.second];
                auto predict_pos = (*a)(key) * this->buffer_ratio;  // normalize to buffer size;
                bufferpiece->model_insert(key,payload,predict_pos);
//                if (mflags)
//                    interchain->remove(interchain->head->next->key);
////                        cout << "i need You, my Lord" <<endl;}
                interchain->put(bufferpiece->key_slots_[0],bufferpiece);
            }
            else{
                // 执行插入
                auto predict_pos = (*a)(key) * this->buffer_ratio;  // normalize to buffer size;
                auto lrukey = bufferpiece->key_slots_[0];
                auto mflags = bufferpiece->model_insert(key,payload,predict_pos);
                if (mflags < 0){
                    interchain->remove( lrukey);
//                        cout << "i need You, my Lord" <<endl;
                }
                interchain->put(bufferpiece->key_slots_[0],bufferpiece);
            }
            /*  // after insert, judge whether to retrain the leaf
            auto lrukey = bufferpiece->key_slots_[0];
            auto predict_pos = (*a)(key) * this->buffer_ratio;  // normalize to buffer size;


            auto mflags = bufferpiece->model_insert(key,payload,predict_pos);
            if (mflags.first)
            {  // findleaf rebuild
                if (mflags.second){
                    interchain->remove( lrukey);
//                        cout << "i need You, my Lord" <<endl;
                }
                a->mlearnedMergeSort(ErrorRecursive, interchain);
                // retrain model, get the new slope and intercept
                int newmodels = internal::make_segmentation( this,a,Error,search_res.first->slotkey,interchain,search_res.second) - link_num;
                this->newmodels += newmodels;
                this->retrain ++;
            }
            else{
                if (mflags.second){
                    interchain->remove( lrukey);
                }
//                auto newlrukey = bufferpiece->key_slots_[0];
//                    if (key != newlrukey)
//                        COUT_THIS("please come, my lovely Lord, key != newlrukey")
                interchain->put(bufferpiece->key_slots_[0],bufferpiece);
            }
            */
        }



        /**
          * statistics information
          **/

        inline size_t leafpieces_count() const{
            if (leaflevel.leafpieces.empty())
                return 0;
            else{
                size_t leaves = leaflevel.leafpieces.size();
                return leaves;
            }
        }

        inline vector<int > each_inner_count() const{
            vector<int > inners(innerlevels.size());
            for (int i = innerlevels.size()-1; i >=0; --i)
                inners[i] = innerlevels[i]->innerpieces.size();
            return inners;
        }

        inline int inner_count() const{
            int innersnum = 0;
            for (int i = innerlevels.size()-1; i >=0; --i)
                innersnum += innerlevels[i]->innerpieces.size();
            return innersnum;
        }

        size_t height() const {
            return innerlevels.size()+1;
        }

        inline void verify_film(film_stats *fstats){
            fstats->leaves = leafpieces_count();
            fstats->innernum = inner_count();
            auto a = each_inner_count();
            fstats->inners.assign( a.begin(), a.end());
            fstats->numlevel = height();
        }


        const map<string, double> show_verify()
        {
            std::map<std::string,double> treeinfo;
            this->verify_film(&vstats);
            this->compute_fill_factor();
            double innodesize = sizeof(key_type) + 8+8;  // startkey, double(slope,intercept)
            double leafnodesize = 2 * sizeof(key_type) + 8 + 8 ; // startkey, endkey, double(slope,intercept),  slotkey(算入了data usage 中，因为在introchain 中 只存储了value)
            double buffernodesize = 12;  // without start key, has: capacity(4), num in buf(4), bitmapsize（4）, (locbitmap, slot key, slot data: 预分配空间)
//            double hashlrusize = sizeof(key_type) + 8 + 16 + 16;   // hash (key,v); doubly-link chain（k，,v, prev,next）
            double hashlrusize = 24 +28;   // hash (key,v, next); doubly-link chain（k(unsigned int)，,v, prev,next）  这个更准确
            double no_hashsize = 16+4;  //prev,next, slot-short int
            double dataV = (valuesize+1)*8;
            double datausage = double(inkeynum*(dataV))/1024/1024;
            double addusage = double((exkeynum+inkeynum)*(1+8)+exkeynum*sizeof(key_type) )/double(1048576);  // exkey, flag, pageID,offset   (pageid-int-8 offset short int), flag for data in memory
            double buffusage = (buffernodesize * vstats.leaves)/double(1048576);
            // addusage: (exkeynum+inkeynum)*(1+8) —— bitmap + pointer;  exkeynum*8 —— pageid(ungined int),offset(unsigned int); exkeynum*sizeof(key_type) ——exkey;
            double indexusge = ((vstats.leaves+this->newmodels)*leafnodesize+vstats.innernum*innodesize + 8 )/double(1048576)+buffusage;  //  leaf nodes, inner nodes, root
            double preusage = reserve_buff*(17)/double(1048576);
            double usepreusage = key_in_buff*(17)/double(1048576);
            double lruusage = double(no_hashsize*(inkeynum-key_in_buff) )/double(1048576);// buffer 不采用local chain, 这个是local lru usage
            // num_of_key_in_model = key in buf.遍历buffe，统计key的数量，类比 compute fill factor
//            double lruusage = double(no_hashsize*inkeynum + hashlrusize * vstats.leaves)/1024/1024;  // buffer 采用 local chain
            treeinfo.insert(pair<std::string,int>("leaves",vstats.leaves));
            treeinfo.insert(pair<std::string,int>("inners",vstats.innernum));
            treeinfo.insert(pair<std::string,int>("levels",vstats.numlevel));
            treeinfo.insert(pair<std::string,double>("datausage",datausage));
            treeinfo.insert(pair<std::string,double>("addusage",addusage));  //key,flag,pageID,offset
            treeinfo.insert(pair<std::string,double>("indexusage",indexusge));
            treeinfo.insert(pair<std::string,double>("lruusage",lruusage));
            treeinfo.insert(pair<std::string,double>("preusage",preusage));
            treeinfo.insert(pair<std::string,double>("usepreusage",usepreusage));
            return treeinfo;
        }

        vector<int> traverse_leaves(){  // traverse the leaf level, find the min leaf and the max leaf
            leaflevel;
            int minslot = leaflevel.leafpieces[0]->slotkey.size();
            int maxslot = leaflevel.leafpieces[0]->slotkey.size();
//            for(auto& v : leaflevel){
//                if (v->slotkey.size() < minslot)
//                    minslot = v->slotkey.size();
//                if (v->slotkey.size() > maxslot)
//                    maxslot = v->slotkey.size();
//            }
            for (unsigned int i = 1; i<leaflevel.leafpieces.size()-1;i++){
                if (leaflevel.leafpieces[i]->slotkey.size() < minslot)
                    minslot = leaflevel.leafpieces[i]->slotkey.size();
                if (leaflevel.leafpieces[i]->slotkey.size() > maxslot)
                    maxslot = leaflevel.leafpieces[i]->slotkey.size();
            }
            int lastslot = leaflevel.leafpieces[leaflevel.leafpieces.size()-1]->slotkey.size();
            cout << " the last slot: " << lastslot << endl;
            vector<int> leafinfo;
            leafinfo.push_back(minslot);
            leafinfo.push_back(maxslot);
            leafinfo.push_back(lastslot);
            return leafinfo;
        }


    };

#pragma pack(push, 1)
    /**
     * 开始定义 FILM 中用的的struct
     */



#pragma pack(pop)

//    typedef filminsert::FILMinsert< key_type, key_type* > filmadatype;
    typedef filminsert::FILMinsert<key_type, key_type*, AlexCompare, std::allocator<std::pair<key_type, key_type*>>, true> filmadatype;

    typedef pair<filmadatype::Leafpiece*,unsigned short int>  filmadalrupair;  //lru 中value 的type ，是一个pair，first is 所属的leaf，second 是slot in leaf
    typedef  adalru::hashLRU <key_type,filmadatype::Leafpiece* , adalru::Node<key_type ,filmadatype::Leafpiece*>* > filmadalrutype;
    typedef adalru::localLRU <unsigned short,key_type* > locallrutype;
    typedef filmstorage::filmdisk<key_type> filmadadisk;
    typedef filmstorage::filmmemory<key_type,key_type*,filmadatype, filmadalrutype,filmadadisk> filmadamemory;
    typedef filmadamemory::memoryusage memoryusage;
    //    adalru::hashLRU <key_type,Leafpiece* , adalru::Node<key_type ,Leafpiece*>* >


    double runtimeevictkeytopage(filmadamemory *filmmem,double totalusemem, filmadatype *filmindex,filmadadisk *diskpage , access_stats *r_stats){
        struct timeval wdt1,wdt2,lt1,lt2;
        double wdtimeuse =0.0;  // 写磁盘的时间
        double ltimeuse = 0.0;

        if (filmmem->evictPoss.size() != 0){
            gettimeofday(&wdt1, NULL);
            for (unsigned int i = 0; i< filmmem->evictPoss.size(); i++){
                auto writeevict =  filmmem->evictPoss[i];
                // evict data, get the tail node, remove the tail of intrachain
                gettimeofday(&lt1, NULL);

                auto evictleaf = filmmem->lru->get_tail1();

                if (evictleaf->value->buffer_piece == NULL){
                    filmadatype::sortPiece * evictbuffer = (filmadatype::sortPiece*) (evictleaf->value);

                    filmmem->lru->remove(evictleaf->key);
                    // the tail of the accessed leaf
                    gettimeofday(&lt2, NULL);
                    ltimeuse += (lt2.tv_sec - lt1.tv_sec) + (double) (lt2.tv_usec - lt1.tv_usec) / 1000000.0;

                    for (int i = 0; i < evictbuffer->slotkey.size();i ++ ){
                        // reduced_evicttable

                        auto writeevict = filmmem->reduced_evictkeytoinpage(evictbuffer->slotkey[i], (key_type*)evictbuffer->slotdata[i],
                                                                            diskpage);

                        evictbuffer->locbitmap[i] = false;
                        evictbuffer->slotdata[i] = (void*)writeevict;

                    }

                    // reduced_evicttable
                }
                else{

                    auto evictslotV = evictleaf->value->intrachain.poptail();   // the tail of the accessed leaf

                    gettimeofday(&lt2, NULL);
                    ltimeuse += (lt2.tv_sec - lt1.tv_sec) + (double) (lt2.tv_usec - lt1.tv_usec) / 1000000.0;
                    auto writeevict = filmmem->reduced_evictkeytoinpage(evictleaf->value->slotkey[evictslotV->key], evictslotV->value,
                                                                        diskpage);
                    evictleaf->value->locbitmap[evictslotV->key] = false;
                    evictleaf->value->slotdata[evictslotV->key] = (void*)writeevict;

                    delete evictslotV;
                    evictslotV = NULL;

                }

            }
            filmindex->inkeynum -= filmmem->evictPoss.size();
            filmindex->exkeynum += filmmem->evictPoss.size();
            filmmem->evictPoss.clear();
            // 批量将inmempages写出磁盘
            int pagenum = filmmem->runtimeevictpagestodisk(diskpage);
            r_stats->writenum += pagenum;
            gettimeofday(&wdt2, NULL);
            wdtimeuse += (wdt2.tv_sec - wdt1.tv_sec) + (double) (wdt2.tv_usec - wdt1.tv_usec) / 1000000.0;
            wdtimeuse -= ltimeuse;
        }
        else{  // just evict from the original table
            double ratio = (totalusemem-filmmem->threshold)/filmmem->threshold  ;
            unsigned int batch_evict = ceil(filmindex->inkeynum*ratio*0.6) + 300 ;
            if (ratio>0.2){
                batch_evict = ceil(filmindex->inkeynum*ratio*0.5);
            }
            gettimeofday(&wdt1, NULL);
            for (unsigned int i = 0; i< batch_evict; i++)   // the number of records to be evicted in batches
            {
                // evict data, get the tail node, remove the tail of intrachain
                gettimeofday(&lt1, NULL);
                auto evictleaf = filmmem->lru->get_tail1();

                if (evictleaf->value->buffer_piece == NULL){
                    filmadatype::sortPiece * evictbuffer = (filmadatype::sortPiece*) (evictleaf->value);

                    filmmem->lru->remove(evictleaf->key);

                    auto evictslotV = evictbuffer->intrachain.poptail();
                    // the tail of the accessed leaf
                    gettimeofday(&lt2, NULL);
                    ltimeuse += (lt2.tv_sec - lt1.tv_sec) + (double) (lt2.tv_usec - lt1.tv_usec) / 1000000.0;




                    for (int i = 0; i < evictbuffer->slotkey.size();i ++ ){
                        // reduced_evicttable

                        auto writeevict = filmmem->reduced_evictkeytoinpage(evictbuffer->slotkey[i], (key_type*)evictbuffer->slotdata[i],
                                                                            diskpage);

                        evictbuffer->locbitmap[i] = false;
                        evictbuffer->slotdata[i] = (void*)writeevict;

                    }

                }
                else{


                    auto evictslotV = evictleaf->value->intrachain.poptail();
                    // the tail of the accessed leaf
                    gettimeofday(&lt2, NULL);

                    auto writeevict = filmmem->reduced_evictkeytoinpage(evictleaf->value->slotkey[evictslotV->key], evictslotV->value,
                                                                        diskpage);
//                    auto writeevict = filmmem->evictkeytoinpage(evictleaf->value->slotkey[evictslotV->key], evictslotV->value,
//                                                                diskpage);

                    evictleaf->value->locbitmap[evictslotV->key] = false;
                    evictleaf->value->slotdata[evictslotV->key] = (void*)writeevict;
                    delete evictslotV;
                    evictslotV = NULL;
                }

            }
            filmindex->inkeynum -= batch_evict;
            filmindex->exkeynum += batch_evict;

            int pagenum = filmmem->runtimeevictpagestodisk(diskpage);
            r_stats->writenum += pagenum;
            gettimeofday(&wdt2, NULL);
            wdtimeuse += (wdt2.tv_sec - wdt1.tv_sec) + (double) (wdt2.tv_usec - wdt1.tv_usec) / 1000000.0;
            wdtimeuse -= ltimeuse;
        }
        r_stats->wdisktime += wdtimeuse;
        r_stats->lrutime += ltimeuse;
        return wdtimeuse;
    }


    double runtimeevictkeytopage2(filmadamemory *filmmem,double totalusemem, filmadatype *filmindex,filmadadisk *diskpage , access_stats *r_stats){

        double wdtimeuse =0.0;  // 写磁盘的时间
//        double ltimeuse = 0.0;

        if (filmmem->reused_index_pageid.size() != 0){
//            gettimeofday(&wdt1, NULL);

            auto reused_pagenum = filmmem->reused_index_pageid.size();   // 判断是否有 reused page
            // 根据 内存占用，判断要写出的数据量，以及要使用的reused_page
            double ratio = (totalusemem-filmmem->threshold)/filmmem->threshold;
            unsigned int batch_page_num = std::min<int>(ceil(filmindex->inkeynum*ratio*0.7)/diskpage->recordnum + 50,reused_pagenum);
            auto wstart_time = std::chrono::high_resolution_clock::now();
            for (unsigned int i = 0; i< batch_page_num; i++){
                // 读取reused page，并且一次写满该page
                unsigned int reused_pageid = filmmem->reused_pageids[i];
                unsigned long reused_index_pageid = filmmem->reused_index_pageid[i];
                // 为该page 创建一个inpage in memory
                auto reused_page = filmmem->createinmempage(reused_pageid,diskpage->pagesize,diskpage->recordnum);
                while (reused_page->freespace){
                    // evict data, get the tail node, remove the tail of intrachain
//                    gettimeofday(&lt1, NULL);
                    auto lstart_time = std::chrono::high_resolution_clock::now();
                    auto evictleaf = filmmem->lru->get_tail1();

                    if (evictleaf->value->buffer_piece == NULL){
                        filmadatype::sortPiece * evictbuffer = (filmadatype::sortPiece*) (evictleaf->value);

//                        ltimeuse += std::chrono::duration_cast<std::chrono::nanoseconds>(
//                                std::chrono::high_resolution_clock::now() - lstart_time)
//                                .count();
                        // 什么时候remove    filmmem->lru->remove(evictleaf->key);
                        /*
                        for (int j = 0; j < evictbuffer->slotkey.size();j ++ ){
                           if (evictbuffer->locbitmap[j]){
                               auto writeevict = filmmem->reduced_evictkeytoinpage(evictbuffer->slotkey[j], (key_type*)evictbuffer->slotdata[j],
                                                                                   diskpage,reused_page, reused_index_pageid);
                               evictbuffer->locbitmap[j] = false;
                               evictbuffer->slotdata[j] = (void*)writeevict;

                               if (!(reused_page->freespace))  //如果还能容纳一个record，那么就向该页写出，如果不能，那么就将当前的inpage写出disk，创建新的inpage
                               {
                                   filmmem->inmempages.emplace_back(reused_page);
                                   if (++i == reused_pagenum)
                                       break;
                                   reused_pageid = filmmem->reused_pageids[i];
                                   reused_index_pageid = filmmem->reused_index_pageid[i];
                                   // 为该page 创建一个inpage in memory
                                   reused_page = filmmem->createinmempage(reused_pageid,diskpage->pagesize,diskpage->recordnum);
                               }
                           }
                        }
                        */

                        // check_exists() and check_whether in memory;
                        for (int k = 0; k < evictbuffer->data_capacity_;k++){
                            if (evictbuffer->check_exists(k)){
                                if (evictbuffer->ALEX_DATA_NODE_FLAG_AT(k)){
                                auto writeevict = filmmem->reduced_evictkeytoinpage(evictbuffer->ALEX_DATA_NODE_KEY_AT(k),
                                                                                    (key_type*)evictbuffer->ALEX_DATA_NODE_PAYLOAD_AT(k),
                                                                                    diskpage,reused_page, reused_index_pageid);
                                evictbuffer->ALEX_DATA_NODE_FLAG_AT(k) = false;
                                evictbuffer->ALEX_DATA_NODE_PAYLOAD_AT(k) = (void*)writeevict;

                                if (!(reused_page->freespace))  //如果还能容纳一个record，那么就向该页写出，如果不能，那么就将当前的inpage写出disk，创建新的inpage
                                {
                                    filmmem->inmempages.emplace_back(reused_page);
                                    if (++i == batch_page_num)
                                        break;
                                    reused_pageid = filmmem->reused_pageids[i];
                                    reused_index_pageid = filmmem->reused_index_pageid[i];
                                    // 为该page 创建一个inpage in memory
                                    reused_page = filmmem->createinmempage(reused_pageid,diskpage->pagesize,diskpage->recordnum);
                                }
                                }
                            }
                        }


                        filmmem->lru->remove(evictleaf->key);
                    }
                    else{
                        auto evictslotV = evictleaf->value->intrachain.poptail();   // the tail of the accessed leaf
//                        ltimeuse += std::chrono::duration_cast<std::chrono::nanoseconds>(
//                                std::chrono::high_resolution_clock::now() - lstart_time)
//                                .count();   // ns
                        auto writeevict = filmmem->reduced_evictkeytoinpage(evictleaf->value->slotkey[evictslotV->key], (key_type*)evictslotV->value,
                                                                            diskpage,reused_page, reused_index_pageid);

                        evictleaf->value->locbitmap[evictslotV->key] = false;
                        evictleaf->value->slotdata[evictslotV->key] = (void*)writeevict;
                        delete evictslotV;
                        evictslotV = NULL;
                    }
                }

                filmmem->inmempages.emplace_back(reused_page);

//                delete reused_page;
            }
            if (filmmem->inmempages.size()>batch_page_num)
                filmmem->inmempages.pop_back();

            filmindex->inkeynum -= (diskpage->recordnum * batch_page_num);
            filmindex->exkeynum += (diskpage->recordnum * batch_page_num);

            if (batch_page_num == reused_pagenum){ //不能再释放所有的reused_pageids,需要判断
                vector<unsigned long >().swap(filmmem->reused_index_pageid);
                vector<pageid_type>().swap(filmmem->reused_pageids);
            }
            else{
                filmmem->reused_index_pageid.assign(filmmem->reused_index_pageid.begin()+batch_page_num,filmmem->reused_index_pageid.end());
                filmmem->reused_pageids.assign(filmmem->reused_pageids.begin()+batch_page_num,filmmem->reused_pageids.end());
            }

            // 需要单独将每个page 写出去磁盘，因为reused 的page 不一定连续
            int pagenum = filmmem->runtimeevictreusedpagestodisk(diskpage);
            r_stats->writenum += pagenum;
            wdtimeuse += std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::high_resolution_clock::now() - wstart_time)
                    .count();   // ns

//            gettimeofday(&wdt2, NULL);
//            wdtimeuse += (wdt2.tv_sec - wdt1.tv_sec) + (double) (wdt2.tv_usec - wdt1.tv_usec) / 1000000.0;
//            wdtimeuse -= ltimeuse;
        }
        else{  // just evict from the original table
            double ratio = (totalusemem-filmmem->threshold)/filmmem->threshold  ;
            unsigned int batch_evict = ceil(filmindex->inkeynum*ratio*0.6) + 300 ;
            if (ratio>0.2){
                batch_evict = ceil(filmindex->inkeynum*ratio*0.5);
            }
//            gettimeofday(&wdt1, NULL);
            auto wstart_time = std::chrono::high_resolution_clock::now();
            for (unsigned int i = 0; i< batch_evict; i++)   // the number of records to be evicted in batches
            {
                // evict data, get the tail node, remove the tail of intrachain
//                gettimeofday(&lt1, NULL);
                auto lstart_time = std::chrono::high_resolution_clock::now();
                auto evictleaf = filmmem->lru->get_tail1();
//                if (evictleaf->key == 309575466846)
//                    cout << "Jesus, i need You" <<endl;
                if (evictleaf->value->buffer_piece == NULL){

                    filmadatype::sortPiece * evictbuffer = (filmadatype::sortPiece*) (evictleaf->value);

//                    ltimeuse += std::chrono::duration_cast<std::chrono::nanoseconds>(
//                            std::chrono::high_resolution_clock::now() - lstart_time)
//                            .count();
                    /*
                    for (int j = 0; j < evictbuffer->slotkey.size();j ++ ){
                        // reduced_evicttable
                        if (evictbuffer->locbitmap[j]){
                            auto writeevict = filmmem->reduced_evictkeytoinpage(evictbuffer->slotkey[j], (key_type*)evictbuffer->slotdata[j],
                                                                                diskpage);
                            evictbuffer->locbitmap[j] = false;
                            evictbuffer->slotdata[j] = (void*)writeevict;
                        }
                    }
                    */

                    for (int k = 0;k < evictbuffer->data_capacity_; k++){
                        if (evictbuffer->check_exists(k) && evictbuffer->ALEX_DATA_NODE_FLAG_AT(k)){
                            auto writeevict = filmmem->reduced_evictkeytoinpage(evictbuffer->ALEX_DATA_NODE_KEY_AT(k),
                                                                                (key_type*)evictbuffer->ALEX_DATA_NODE_PAYLOAD_AT(k),
                                                                                diskpage);
                            evictbuffer->ALEX_DATA_NODE_FLAG_AT(k) = false;
                            evictbuffer->ALEX_DATA_NODE_PAYLOAD_AT(k) = (void*)writeevict;

                        }
                    }

                    filmmem->lru->remove(evictleaf->key);
                }
                else{
//                    gettimeofday(&lt1, NULL);
                    auto evictslotV = evictleaf->value->intrachain.poptail();
                    // the tail of the accessed leaf
//                    ltimeuse += std::chrono::duration_cast<std::chrono::nanoseconds>(
//                            std::chrono::high_resolution_clock::now() - lstart_time)
//                            .count();   // ns

                    auto writeevict = filmmem->reduced_evictkeytoinpage(evictleaf->value->slotkey[evictslotV->key], evictslotV->value,
                                                                        diskpage);

                    evictleaf->value->locbitmap[evictslotV->key] = false;
                    evictleaf->value->slotdata[evictslotV->key] = (void*)writeevict;
                    delete evictslotV;
                    evictslotV = NULL;
                }

            }
            filmindex->inkeynum -= batch_evict;
            filmindex->exkeynum += batch_evict;

            int pagenum = filmmem->runtimeevictpagestodisk(diskpage);
            r_stats->writenum += pagenum;
//            gettimeofday(&wdt2, NULL);
//            wdtimeuse += (wdt2.tv_sec - wdt1.tv_sec) + (double) (wdt2.tv_usec - wdt1.tv_usec) / 1000000.0;
            wdtimeuse += std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::high_resolution_clock::now() - wstart_time)
                    .count();
//            wdtimeuse -= ltimeuse;
        }

        wdtimeuse = wdtimeuse/1e9;
//        ltimeuse = ltimeuse/1e9;
        r_stats->wdisktime += wdtimeuse;
//        r_stats->lrutime += ltimeuse;
        return wdtimeuse;
    }



    // 采用unordered_map 实现 dict 的功能
    template<class v_type>   // int bplustree:  pair< oribttype::leaf_node*,pair<short int,int>
    class prepass_dict{
    public:
        int pagenum = 0;
//        std::unordered_map<int,vector<v_type>> prepass;
        std::map< pageid_type,vector<v_type>> prepass;
        struct dict_find {
            /// find result flags
            bool flags;
            vector<v_type> *findinfo;

            inline dict_find()
                    : flags() {}
        };

        void insert_dict(pageid_type pageid,v_type keyinfo){  //keyinfo(first: offset in page, evicted posi pair<int,short>*)
//            typename std::unordered_map<int,vector<v_type>>::const_iterator got = prepass.find(pageid);
            typename std::map<pageid_type,vector<v_type>>::const_iterator got = prepass.find(pageid);
            if (got == prepass.end())
            {
                prepass[pageid].push_back(keyinfo);
                pagenum +=1;
            }
            else
                prepass[pageid].push_back(keyinfo);
        }

        void insert_dict(pageid_type pageid,pageOff_type off, pair<pageid_type ,pageOff_type> * evictposi){  //keyinfo(first: offset in page, evicted posi pair<int,short>*)
//            typename std::unordered_map<int,vector<v_type>>::const_iterator got = prepass.find(pageid);
            typename std::map<pageid_type,vector<v_type>>::const_iterator got = prepass.find(pageid);
            if (got == prepass.end())
            {
                prepass[pageid].push_back(off);
                pagenum +=1;
            }
            else
                prepass[pageid].push_back(off);
        }


        dict_find find_in_dict(pageid_type pageid){
            dict_find res;
//            typename std::unordered_map<int,vector<v_type>>::const_iterator got = prepass.find(pageid);
            typename std::map<pageid_type,vector<v_type>>::const_iterator got = prepass.find(pageid);
            if (got == prepass.end()){
                std::cout<<"not find the key in prepass, my lovely Lord, waiting for you! "<<endl;
                res.flags = false;
                return res;
            }
            else
            {
                res.flags = true;
                res.findinfo = got;
                return res;
            }
        }
    };

    typedef prepass_dict< pair< filmadatype::Leafpiece* ,pair< lruOff_type, pageOff_type > > > pre_dict ;

// in the range query, judge whether the keys in each overlapped leaf in memory or in disk
    void range_is_in_or_out(int start, int end, filmadatype::Leafpiece* curleaf ,vector<pair<key_type,key_type*>> *totalres , pre_dict* dict,int sizevalue ,filmadamemory *filmmem,access_stats *r_stats){
        struct timeval lt1, lt2;
        double ltimeuse;
        for (int k = start;k < end;k++)
        {
            bool flag = curleaf->locbitmap[k];
            if (flag){
                pair<key_type,key_type*> resitem;
                resitem.first = curleaf->slotkey[k];
                auto finddata = (adalru::Node<lruOff_type, key_type* >*) curleaf->slotdata[k];
                resitem.second = finddata->value;
                totalres->emplace_back(resitem);
                gettimeofday(&lt1, NULL);

                curleaf->intrachain.moveTohead(finddata);// update intrachain
                filmmem->lru->put(curleaf->startkey,curleaf);                 // 更新global lru
                gettimeofday(&lt2, NULL);
                ltimeuse = (lt2.tv_sec - lt1.tv_sec) + (double) (lt2.tv_usec - lt1.tv_usec) / 1000000.0;
                r_stats->lrutime += ltimeuse;
            }
            else{  // prepass into dict, that access disk
                auto evictpos = (pair<pageid_type ,pageOff_type>*)curleaf->slotdata[k];
                pair<pageid_type ,pageOff_type>* diskpos = evictpos;
                pageid_type pageid = diskpos->first;
                auto off = diskpos->second;
                //dict 中插入的是 pageid，所属leafpiece，（在leaf 中的slot 位置，保存evictposi的地址（ pointer ）——用于后期存储被移出key 的外存地址）
                dict->insert_dict(pageid, pair< filmadatype::Leafpiece*,pair<lruOff_type,pageOff_type > >(curleaf,pair<lruOff_type,pageOff_type>(k,off)));

//                cout<<"my Lord, my heart is designed for You!"<<endl;
            }
        }
    }

    void htap_range_is_in_or_out(int start, int end, filmadatype::Leafpiece* curleaf ,vector<pair<key_type,key_type*>> *totalres , pre_dict* dict,
                                 int sizevalue ,filmadamemory *filmmem,access_stats *r_stats,int reduce_factor){
        struct timeval lt1, lt2;
        double ltimeuse;
        for (int k = start;k < end;k++)
        {
            bool flag = curleaf->locbitmap[k];
            if (flag){
                pair<key_type,key_type*> resitem;
                resitem.first = curleaf->slotkey[k];
                auto finddata = (adalru::Node<lruOff_type, key_type* >*) curleaf->slotdata[k];
                resitem.second = finddata->value;
                totalres->emplace_back(resitem);
                gettimeofday(&lt1, NULL);

                curleaf->intrachain.moveTohead(finddata);// update intrachain
                filmmem->lru->put(curleaf->startkey,curleaf);                 // 更新global lru
                gettimeofday(&lt2, NULL);
                ltimeuse = (lt2.tv_sec - lt1.tv_sec) + (double) (lt2.tv_usec - lt1.tv_usec) / 1000000.0;
                r_stats->lrutime += ltimeuse;
            }
            else{  // prepass into dict, that access disk

                pair<pageid_type, pageOff_type>  reduced_writeevict (0,0);

                auto short_reduced_evict = (unsigned long) curleaf->slotdata[k];  //writeevict 指向的是 pageid and offset
                auto int_reduced_evict = short_reduced_evict/2.0;  //308765
                auto fraction = modf(int_reduced_evict,&int_reduced_evict);
                auto int_high_low = filmmem->reduced_evicttable[int_reduced_evict];
                // 根据reduced_evict 算出（读出） page id，算出读出offset
                if (fraction){
                    reduced_writeevict.second = int_high_low&0xFFFF;
//                                unsigned short hi = hhh>>16 ;// 高16位
//                                unsigned short lo = hhh&0xFFFF ;// 低16位
                }
                else
                    reduced_writeevict.second =  int_high_low>>16 ;   // 高16位
                unsigned long index_pageid = floor(int_reduced_evict/reduce_factor)*reduce_factor;
                reduced_writeevict.first = filmmem->reduced_evicttable[index_pageid];
                // 将 querykey merge到memory，修改 bitmap，修改num in page
//                            auto numinpage = memoryfilm.reduced_evicttable[index_pageid+1];
//                filmmem->reduced_evicttable[index_pageid+1] -= 1;


                pageid_type pageid = reduced_writeevict.first;
                auto off = reduced_writeevict.second;
                //dict 中插入的是 pageid，所属leafpiece，（在leaf 中的slot 位置，保存evictposi的地址（ pointer ）——用于后期存储被移出key 的外存地址）
                dict->insert_dict(pageid, pair< filmadatype::Leafpiece*,pair<lruOff_type,pageOff_type > >(curleaf,pair<lruOff_type,pageOff_type>(k,off)));

//                cout<<"my Lord, my heart is designed for You!"<<endl;
            }
        }
    }

    void htap_range_is_in_or_out(int start, int end, filmadatype::sortPiece* curleaf ,vector<pair<key_type,key_type*>> *totalres , pre_dict* dict,int sizevalue ,filmadamemory *filmmem,access_stats *r_stats,int reduce_factor){
        struct timeval lt1, lt2;
        double ltimeuse;
        for (int k = start;k < end;k++)
        {
            if (curleaf->check_exists(k))  // 如果不是gap
            {
                bool flag = curleaf->ALEX_DATA_NODE_FLAG_AT(k);
                if (flag){
                    pair<key_type,key_type*> resitem;
                    resitem.first = curleaf->ALEX_DATA_NODE_KEY_AT(k);
                    auto finddata = (key_type*) curleaf->ALEX_DATA_NODE_PAYLOAD_AT(k);
                    resitem.second = finddata;
                    totalres->emplace_back(resitem);
                    gettimeofday(&lt1, NULL);

                    filmmem->lru->put(curleaf->ALEX_DATA_NODE_KEY_AT(0),curleaf);                 // 更新global lru
                    gettimeofday(&lt2, NULL);
                    ltimeuse = (lt2.tv_sec - lt1.tv_sec) + (double) (lt2.tv_usec - lt1.tv_usec) / 1000000.0;
                    r_stats->lrutime += ltimeuse;
                }
                else{  // prepass into dict, that access disk
//                    auto evictpos = (pair<pageid_type ,pageOff_type>*)curleaf->ALEX_DATA_NODE_PAYLOAD_AT(k);

                    pair<pageid_type, pageOff_type>  reduced_writeevict (0,0);

                    auto short_reduced_evict = (unsigned long) curleaf->ALEX_DATA_NODE_PAYLOAD_AT(k);  //writeevict 指向的是 pageid and offset
                    auto int_reduced_evict = short_reduced_evict/2.0;
                    auto fraction = modf(int_reduced_evict,&int_reduced_evict);
                    auto int_high_low = filmmem->reduced_evicttable[int_reduced_evict];
                    // 根据reduced_evict 算出（读出） page id，算出读出offset
                    if (fraction){
                        reduced_writeevict.second = int_high_low&0xFFFF;
//                                unsigned short hi = hhh>>16 ;// 高16位
//                                unsigned short lo = hhh&0xFFFF ;// 低16位
                    }
                    else
                        reduced_writeevict.second =  int_high_low>>16 ;   // 高16位
                    unsigned long index_pageid = floor(int_reduced_evict/reduce_factor)*reduce_factor;
                    reduced_writeevict.first = filmmem->reduced_evicttable[index_pageid];
                    // 将 querykey merge到memory，修改 bitmap，修改num in page
//                            auto numinpage = memoryfilm.reduced_evicttable[index_pageid+1];
//                filmmem->reduced_evicttable[index_pageid+1] -= 1;


                    pageid_type pageid = reduced_writeevict.first;
                    auto off = reduced_writeevict.second;
                    //dict 中插入的是 pageid，所属leafpiece，（在leaf 中的slot 位置，保存evictposi的地址（ pointer ）——用于后期存储被移出key 的外存地址）
                    dict->insert_dict(pageid, pair< filmadatype::Leafpiece*,pair<lruOff_type,pageOff_type > >(curleaf,pair<lruOff_type,pageOff_type>(k,off)));

//                cout<<"my Lord, my heart is designed for You!"<<endl;
                }
            }

        }
    }


    int range_prepass(filmadamemory *filmmem,filmadatype::range_res_find *index_result,vector<pair<key_type,key_type*>> *result,pre_dict* dict,filmadatype *filmindex,access_stats *r_stats){  // 将range query 由btree 得到的信息，进行怕热pass，遍历find leaf中的每一个leaf，如果slotdata.first = true, 则直接放入reslut 中，否则将其放入prepassdict 中

        int leafnum = index_result->findleafs.size();
        if (leafnum == 1){  //all the keys in range belonging to a certain leaf
            range_is_in_or_out(index_result->firstslot,index_result->lastslot+1,index_result->findleafs[0],result,dict,filmindex->valuesize,filmmem,r_stats);
        }
        else {  // the first leaf and last leaf need to be deal with seperately
            // deal with firstleaf
            //index_result->findleafs[0];
            range_is_in_or_out(index_result->firstslot,index_result->findleafs[0]->slotkey.size(),index_result->findleafs[0],result,dict,filmindex->valuesize,filmmem,r_stats);
            // deal with the middle leafs
            for (int i=1; i < leafnum-1; i++){
                range_is_in_or_out(0,index_result->findleafs[i]->slotkey.size(),index_result->findleafs[i],result,dict,filmindex->valuesize,filmmem,r_stats);
            }
            //deal with lastleaf
            //index_result->findleafs[leafnum-1];
            range_is_in_or_out(0,index_result->lastslot+1,index_result->findleafs[leafnum-1],result,dict,filmindex->valuesize,filmmem,r_stats);
        }
        return 0;
    }

    int htap_range_prepass(filmadamemory *filmmem,filmadatype::htaprange_res_find *index_result,vector<pair<key_type,key_type*>> *result,pre_dict* dict,filmadatype *filmindex,access_stats *r_stats, int reduce_factor){  // 将range query 由btree 得到的信息，进行怕热pass，遍历find leaf中的每一个leaf，如果slotdata.first = true, 则直接放入reslut 中，否则将其放入prepassdict 中

        int leafnum = index_result->findleafs.size();
        if (leafnum == 1){  //all the keys in range belonging to a certain leaf
//            cout << "Lord, this is chao" << endl;
            htap_range_is_in_or_out(index_result->slotleaf0,index_result->slotleaf1+1,index_result->findleafs[0],result,dict,filmindex->valuesize,filmmem,r_stats,reduce_factor);
        }
        else {  // the first leaf and last leaf need to be deal with seperately
            // deal with firstleaf
//            cout << "Lord, this is chao" << endl;
            htap_range_is_in_or_out(index_result->slotleaf0,index_result->findleafs[0]->slotkey.size(),index_result->findleafs[0],result,dict,filmindex->valuesize,filmmem,r_stats,reduce_factor);
            // deal with the middle leafs
            for (int i=1; i < leafnum-1; i++){
                htap_range_is_in_or_out(0,index_result->findleafs[i]->slotkey.size(),index_result->findleafs[i],result,dict,filmindex->valuesize,filmmem,r_stats,reduce_factor);
            }
            //deal with lastleaf
//            cout << "Lord, this is chao" << endl;
            htap_range_is_in_or_out(0,index_result->slotleaf1+1,index_result->findleafs[leafnum-1],result,dict,filmindex->valuesize,filmmem,r_stats,reduce_factor);
        }

        int buffernum = index_result->findbuffers.size();
        if (buffernum > 0){
            index_result->slotbuffer0 = std::max<int>(0,index_result->slotbuffer0);
            if (buffernum == 1){  //all the keys in range belonging to a certain leaf
//                cout << "Lord, this is chao" << endl;
//                if (index_result->findbuffers[0]->num_keys_ > 0)
                if (index_result->slotbuffer0!=index_result->slotbuffer1) {
                    // 找到slotbuffer0 开始的第一个 valid pos
                    auto buffer = (filmadatype::sortPiece *) index_result->findbuffers[0];
                    auto start_buffer0 = buffer->find_first_1(index_result->slotbuffer0);  // 需要考虑
//                    htap_range_is_in_or_out(index_result->slotbuffer0, index_result->slotbuffer1,
//                                            buffer, result, dict,
//                                            filmindex->valuesize, filmmem, r_stats, reduce_factor);
                    htap_range_is_in_or_out(start_buffer0, index_result->slotbuffer1,
                                            buffer, result, dict,
                                            filmindex->valuesize, filmmem, r_stats, reduce_factor);
                }
            }
            else {  // the first leaf and last leaf need to be deal with seperately
                // deal with firstleaf
//                cout << "Lord, this is chao" << endl;
                filmadatype::sortPiece* buffer = (filmadatype::sortPiece*) (index_result->findbuffers[0]);
                auto start_buffer0 = buffer->find_first_1(index_result->slotbuffer0); // 需要考虑//  auto start_buffer0 = buffer->find_first_1();   //不需要考虑
//                auto start_buffer00 = buffer->find_first_1_incorrect(index_result->slotbuffer0); //
                auto end_buffer0 = buffer->find_last_1()+1;   // 需要考虑
//                htap_range_is_in_or_out(index_result->slotbuffer0,buffer->data_capacity_,buffer,result,dict,filmindex->valuesize,filmmem,r_stats,reduce_factor);
                htap_range_is_in_or_out(start_buffer0,end_buffer0,(filmadatype::sortPiece*)index_result->findbuffers[0],result,dict,filmindex->valuesize,filmmem,r_stats,reduce_factor);
                // deal with the middle leafs
                for (int i=1; i < buffernum-1; i++){
                    buffer = (filmadatype::sortPiece*) (index_result->findbuffers[i]);
                    // 找到该buffer 的第一个valid record 的位置
                    auto start_buffer = buffer->find_first_1();  //需要考虑
                    auto end_buffer = buffer->find_last_1()+1;  // 需要考虑
//                    htap_range_is_in_or_out(0,buffer->data_capacity_,index_result->findbuffers[i],result,dict,filmindex->valuesize,filmmem,r_stats,reduce_factor);
                    htap_range_is_in_or_out(start_buffer,end_buffer,index_result->findbuffers[i],result,dict,filmindex->valuesize,filmmem,r_stats,reduce_factor);
//                    COUT_THIS("thank You, my lovely Lord");
                }
                //deal with lastleaf
//                cout << "Lord, this is chao" << endl;
                buffer = (filmadatype::sortPiece*)(index_result->findbuffers[buffernum-1]);
                auto start_buffer = buffer->find_first_1();  // 需要考虑。找到该buffer 的第一个valid record 的位置
//                auto end_buffer = buffer->find_last_1();  // 不需要考虑
//                htap_range_is_in_or_out(0,index_result->slotbuffer1,buffer,result,dict,filmindex->valuesize,filmmem,r_stats,reduce_factor);
                htap_range_is_in_or_out(start_buffer,index_result->slotbuffer1,buffer,result,dict,filmindex->valuesize,filmmem,r_stats,reduce_factor);
            }
        }

        return 0;
    }


    // 一次至少读取1024k 的数据
    int range_read_from_disk(vector<pair<key_type,key_type*>> *result,pre_dict* dict,filmadamemory *filmmem,
                             filmadatype *filmindex,filmadadisk *pagedisk, access_stats *r_stats, key_type* payload) {
        // 打开磁盘文件
        struct timeval lt1, lt2, dt1, dt2, rdt1, rdt2;
        double ltimeuse = 0.0;  // 读磁盘的时间，写磁盘的时间
        double dtimeuse, rdtimeuse;
        unsigned long cross = dict->prepass.rbegin()->first - dict->prepass.begin()->first + 1;
//        cout<< cross << " ";
//        r_stats->pagecross += dict->pagenum;
        int fd;
        key_type *buf;

        pair<key_type, key_type *> res;

        unsigned long int fix_buf_size = pagedisk->pagesize * 8;  // 每个磁盘页的大小


        fd = open(pagedisk->pagefile, O_RDWR | O_DIRECT, 0755);
        if (cross * pagedisk->pagesize <= 131072)// 262144 2048k 判断一下，如果最大page 和最小page，相差超过1024k 131072，则采用seek的形式进行读取，
        {
            unsigned long int buf_size = fix_buf_size * cross;  // 一次性读入
            int ret = posix_memalign((void **) &buf, 512, buf_size);
            memset(buf, 'c', buf_size);

            gettimeofday(&rdt1, NULL);
            pageid_type firstp = dict->prepass.begin()->first;
            unsigned long aoffset = firstp * fix_buf_size;
            ret = pread(fd, buf, buf_size, static_cast<unsigned long> (aoffset)); // 从 prepass 的第一个页开始，读取所有跨着的页
            gettimeofday(&rdt2, NULL);
            rdtimeuse = (rdt2.tv_sec - rdt1.tv_sec) + (double) (rdt2.tv_usec - rdt1.tv_usec) / 1000000.0;
            r_stats->rdisktime += rdtimeuse;
            gettimeofday(&dt1, NULL);
            for (auto &v: dict->prepass) {
                // buf 中为所有要读取的数据
                pageid_type pageid = v.first;
                int diff = pageid- dict->prepass.begin()->first;
                unsigned long index_pageid = v.first*pagedisk->reduce_factor;
//                auto iiii = filmmem->reduced_evicttable[index_pageid];
//                if (pageid == 3027)
//                    cout << "Jesus, please teach me!" << endl;
                unsigned long int offset = diff * (pagedisk->recordsize * pagedisk->recordnum);// 把当前页的数据 单独出来
                for (int k = 0; k <
                                v.second.size(); k++) {   //v.second.size()  是要从当前页中读取多少个data， v.second[i] 是每个data 的相关信息是一个pair，v.second[i].first 是所属的leaf， v.second[i].second。first是slot in leaf，second 是 evictposi指向page id 和off
//                    auto slot1 = v.second[k].second.first;
                    auto slot = v.second[k].second.first;
                    pageOff_type off = v.second[k].second.second;   // offset in page
//                    auto writeevict = v.second[k].second.second;


                    if (v.first == filmmem->inpage->pageid) {
                        //vector<key_type> inmemdata = filmmem->inpage->inmemdata;
                        key_type *inmemdata = filmmem->inpage->inmemdata;
                        res.first = inmemdata[off];
                        res.second = new key_type[filmindex->valuesize];
                        for (int i = 0; i < filmindex->valuesize; i++) {
                            res.second[i] = inmemdata[off + 1 + i];
                        }
                        delete [] res.second;

                    } else {
                        res.first = buf[offset + off];
                        res.second = new key_type[filmindex->valuesize];
                        for (int i = 0; i < filmindex->valuesize; i++) {
                            res.second[i] = buf[offset + off + 1 + i];
                        }
                        delete [] res.second;
                    }
                    // 修改该页对应的page info
                    auto recordnum = off/pagedisk->recordsize;  // compute the key's rank in the page
                    filmmem->unset_bit(recordnum, index_pageid);
                    filmmem->reduced_evicttable[index_pageid+1] -= 1;
                    filmindex->inkeynum++;
                    filmindex->exkeynum--;

                    result->push_back(res);

                    if (v.second[k].first->buffer_piece != NULL){
                        filmadatype::Leafpiece *curpiece = v.second[k].first;
                        if (res.first != curpiece->slotkey[slot]) {
                            cout << "i need You, my Lord, Thank You for all the thing!  not find "<< res.first  << endl;
                        }
                        // 将res 插入到curpiece 的intrachain，// 将被访问的record 的数据放入内存中，即其所属的leaf，并修改其slotdata 的flag 为 true
                        gettimeofday(&lt1, NULL);
//                        curpiece->slotdata[slot] = curpiece->intrachain.put(slot, res.second);
                        curpiece->slotdata[slot] = curpiece->intrachain.put(slot, payload);
                        curpiece->locbitmap[slot] = true;

                        // 将curpiece 插入到interchain 中，即mem->lru中
                        filmmem->lru->put(curpiece->startkey, curpiece);// 更新 interchain
                        gettimeofday(&lt2, NULL);
                        ltimeuse += (lt2.tv_sec - lt1.tv_sec) + (double) (lt2.tv_usec - lt1.tv_usec) / 1000000.0;
                    }
                    else{
                        filmadatype::sortPiece *curpiece = (filmadatype::sortPiece *)v.second[k].first;
                        /*if (res.first != curpiece->ALEX_DATA_NODE_KEY_AT(slot)) {
                            cout << "i need You, my Lord, Thank You for all the thing! not find " << res.first << endl;
                        }*/
                        // 将res 插入到curpiece 的intrachain，// 将被访问的record 的数据放入内存中，即其所属的leaf，并修改其slotdata 的flag 为 true
                        gettimeofday(&lt1, NULL);
//                        curpiece->ALEX_DATA_NODE_PAYLOAD_AT(slot) = res.second;
                        curpiece->ALEX_DATA_NODE_PAYLOAD_AT(slot) = payload;
                        curpiece->ALEX_DATA_NODE_FLAG_AT(slot) = true;

                        // 将curpiece 插入到interchain 中，即mem->lru中
                        filmmem->lru->put(curpiece->ALEX_DATA_NODE_KEY_AT(0), curpiece);// 更新 interchain
                        gettimeofday(&lt2, NULL);
                        ltimeuse += (lt2.tv_sec - lt1.tv_sec) + (double) (lt2.tv_usec - lt1.tv_usec) / 1000000.0;
                    }

                }
                // check whether to merge key // 修改该页对应的page info，根据该info 判断是否merge data
                if (filmmem->check_merge(index_pageid) == true){
                    if (pageid != filmmem->inpage->pageid){
//                        auto numinpage = filmmem->reduced_evicttable[index_pageid+1];
                        filmmem->reused_pageids.emplace_back(pageid);
                        filmmem->num_page_merge += 1;
                        filmmem->reused_index_pageid.emplace_back(index_pageid);
                        auto pagenum = filmmem->reduced_evicttable[index_pageid+1];
                        filmmem->page_merge(buf,pagedisk->recordnum, index_pageid,offset);
                        pagenum = filmmem->reduced_evicttable[index_pageid+1];
                        filmindex->exkeynum -= filmmem->reduced_evicttable[index_pageid+1];
                        filmindex->inkeynum += filmmem->reduced_evicttable[index_pageid+1];
                        filmmem->reduced_evicttable[index_pageid+1] = 0;
                    }
//                    else{
//                        // record the index_pageid, that the position in the reduced_evictedtable, which can be reused.
//                        filmmem->reused_pageids.emplace_back(pageid);
//                        filmmem->num_page_merge += 1;
//                        filmmem->reused_index_pageid.emplace_back(index_pageid);
//                        filmmem->page_merge(filmmem->inpage->inmemdata,pagedisk->recordnum, index_pageid);
//                        filmmem->inpage->pageid +=1;
//                        filmmem->inpageid ++;
//                        filmmem->inpage->recordnum = 0;
//                        filmmem->inpage->freespace = pagedisk->pagesize;
//                        filmmem->inpage->usedspace = 0;
//
//                        filmindex->exkeynum -= filmmem->reduced_evicttable[index_pageid+1];
//                        filmindex->inkeynum += filmmem->reduced_evicttable[index_pageid+1];
//                        filmmem->reduced_evicttable[index_pageid+1] = 0;
//                    }
                }

            }
            gettimeofday(&dt2, NULL);
            dtimeuse = (dt2.tv_sec - dt1.tv_sec) + (double) (dt2.tv_usec - dt1.tv_usec) / 1000000.0;
            r_stats->disktime += dtimeuse;
            r_stats->disktime -= ltimeuse;
            r_stats->lrutime += ltimeuse;

        } else {
            int ret = posix_memalign((void **) &buf, 512, fix_buf_size);
            memset(buf, 'c', fix_buf_size);
            // for loop to read the data of each page needed
            pageid_type abslotepageid = 0;
            gettimeofday(&dt1, NULL);
            for (auto &v: dict->prepass) {
                pageid_type pageid = v.first;  // v.first 是 pageid， v.second 是 这个page 要读取的数据的信息，
                unsigned long index_pageid = pageid*pagedisk->reduce_factor;
//                auto iiii = filmmem->reduced_evicttable[index_pageid];
//                if (pageid == 3027)
//                    cout << "Jesus, plese teach me!" << endl;
                if (pageid == filmmem->inpage->pageid) {
//                    cout << " Jesus, Son of David, please have pity on me, read from in-memory page" << endl;
                    key_type *inmemdata = filmmem->inpage->inmemdata;

                    for (int vi = 0; vi <
                                     v.second.size(); vi++) {   //v.second.size()  是要从当前页中读取多少个data， v.second[i] 是每个data 的相关信息是一个pair，v.second[i].first 是所属的leaf， v.second[i].second。first是slot in leaf，second 是 evictposi指向page id 和off
                        auto slot = v.second[vi].second.first;   // v.first 是 page id，v.second 是要从该页中读进来的每个record 的信息，是一个pair，
                        // 该pair 的first 是所属的leaf，该pair 的second 是一个pair，
                        // 包括：first 是在该leaf 中的slot，second 是page 信息，包括page id 和在该page 中的offset
                        short int off = v.second[vi].second.second;  //offset in page
//                        auto writeevict = v.second[vi].second.second;

                        if (v.second[vi].first->buffer_piece != NULL){

                            filmadatype::Leafpiece *curpiece = v.second[vi].first;
                            res.first = inmemdata[off];
                            res.second = new key_type[filmindex->valuesize];
                            for (int ki = 0; ki < filmindex->valuesize; ki++) {
                                res.second[ki] = inmemdata[off + 1 + ki];
                            }// check page merge
                            delete [] res.second;
                            auto recordnum = off/pagedisk->recordsize;  // compute the key's rank in the page
//                        bool former_flag = filmmem->check_exists(recordnum,index_pageid);
//                        if (former_flag == false){
////                auto numinpage = memoryfilm->pageinfo[index_pageid+1];
//                            cout << "i need You, my Lord, the flag is wrong" << endl;
//                        }
                            filmmem->unset_bit(recordnum, index_pageid);
                            filmmem->reduced_evicttable[index_pageid+1] -= 1;
//                        if (former_flag == false){
//                            cout << "i need You, my Lord, the flag is wrong" << endl;
//                        }

                            filmindex->inkeynum++;
                            filmindex->exkeynum--;
                            result->push_back(res);

                            /*if (res.first != curpiece->slotkey[slot]) {
                                cout << "i need You, my Lord, Thank You for all the thing! not find "<<res.first << endl;
                            }*/
                            // 将res 插入到curpiece 的intrachain，// 将被访问的record 的数据放入内存中，即其所属的leaf，并修改其slotdata 的flag 为 true
                            gettimeofday(&lt1, NULL);
//                            curpiece->slotdata[slot] = curpiece->intrachain.put(slot, res.second);
                            curpiece->slotdata[slot] = curpiece->intrachain.put(slot, payload);
                            curpiece->locbitmap[slot] = true;
                            // 将curpiece 插入到interchain 中，即mem->lru中
                            filmmem->lru->put(curpiece->startkey, curpiece);// 更新 interchain
                            gettimeofday(&lt2, NULL);
                            ltimeuse += (lt2.tv_sec - lt1.tv_sec) + (double) (lt2.tv_usec - lt1.tv_usec) / 1000000.0;
//                      cout<< " i need You, my Lord!!! we all need You! "<< endl;
                        }
                        else{

                            filmadatype::sortPiece *curpiece = (filmadatype::sortPiece *)v.second[vi].first;

                            res.first = inmemdata[off];
                            res.second = new key_type[filmindex->valuesize];
                            for (int ki = 0; ki < filmindex->valuesize; ki++) {
                                res.second[ki] = inmemdata[off + 1 + ki];
                            }// check page merge
                            delete [] res.second;
                            auto recordnum = off/pagedisk->recordsize;  // compute the key's rank in the page

                            filmmem->unset_bit(recordnum, index_pageid);
                            filmmem->reduced_evicttable[index_pageid+1] -= 1;
//                        if (former_flag == false){
//                            cout << "i need You, my Lord, the flag is wrong" << endl;
//                        }

                            filmindex->inkeynum++;
                            filmindex->exkeynum--;
                            result->push_back(res);

                            /*if (res.first != curpiece->ALEX_DATA_NODE_KEY_AT(slot)) {
                                cout << "i need You, my Lord, Thank You for all the thing! not find "<<res.first << endl;
                            }*/
                            // 将res 插入到curpiece 的intrachain，// 将被访问的record 的数据放入内存中，即其所属的leaf，并修改其slotdata 的flag 为 true
                            gettimeofday(&lt1, NULL);

//                            curpiece->ALEX_DATA_NODE_PAYLOAD_AT(slot) = res.second;
                            curpiece->ALEX_DATA_NODE_PAYLOAD_AT(slot) = payload;
                            curpiece->ALEX_DATA_NODE_FLAG_AT(slot) = true;

                            // 将curpiece 插入到interchain 中，即mem->lru中
                            filmmem->lru->put(curpiece->ALEX_DATA_NODE_KEY_AT(0), curpiece);// 更新 interchain
                            gettimeofday(&lt2, NULL);
                            ltimeuse += (lt2.tv_sec - lt1.tv_sec) + (double) (lt2.tv_usec - lt1.tv_usec) / 1000000.0;

                        }


                    }

                }
                else {

                    unsigned long seekdis = fix_buf_size * (pageid - abslotepageid);
                    lseek(fd, seekdis, SEEK_CUR);
                    ret = read(fd, buf, fix_buf_size);
                    abslotepageid = pageid + 1;

                    //buf 中为所有属于这个页的数据，根据v.second 包含的offset in page 的信息，将对应的值读出来，并放入leaf 中
                    for (int vi = 0; vi <
                                     v.second.size(); vi++) {   //v.second.size()  是要从当前页中读取多少个data， v.second[i] 是每个data 的相关信息是一个pair，v.second[i].first 是所属的leaf， v.second[i].second。first是slot in leaf，second 是 evictposi指向page id 和off
                        auto slot = v.second[vi].second.first;
                        short int off = v.second[vi].second.second;  //offset in page
//                        auto writeevict = v.second[vi].second.second;
//                        filmadatype::Leafpiece *curpiece = v.second[vi].first;
                        res.first = buf[off];
                        res.second = new key_type[filmindex->valuesize];
                        for (int ki = 0; ki < filmindex->valuesize; ki++) {
                            res.second[ki] = buf[off + 1 + ki];
                        }
                        delete [] res.second;
                        auto recordnum = off/pagedisk->recordsize;  // compute the key's rank in the page
                        filmmem->unset_bit(recordnum, index_pageid);
                        filmmem->reduced_evicttable[index_pageid+1] -= 1;
//           check page merge
                        filmindex->inkeynum++;
                        filmindex->exkeynum--;
                        result->push_back(res);

                        if (v.second[vi].first->buffer_piece != NULL){
                            filmadatype::Leafpiece *curpiece = v.second[vi].first;
                            // 将res 插入到curpiece 的intrachain，// 将被访问的record 的数据放入内存中，即其所属的leaf，并修改其slotdata 的flag 为 true
                            gettimeofday(&lt1, NULL);
//                            curpiece->slotdata[slot] = curpiece->intrachain.put(slot, res.second);
                            curpiece->slotdata[slot] = curpiece->intrachain.put(slot, payload);
                            curpiece->locbitmap[slot] = true;
                            // 将curpiece 插入到interchain 中，即mem->lru中
                            filmmem->lru->put(curpiece->startkey, curpiece);// 更新 interchain
                            gettimeofday(&lt2, NULL);
                            ltimeuse += (lt2.tv_sec - lt1.tv_sec) + (double) (lt2.tv_usec - lt1.tv_usec) / 1000000.0;
//                cout<< " i need You, my Lord!!! we all need You! "<< endl;
                            /*if (res.first != curpiece->slotkey[slot]) {
                                cout << "i need You, my Lord, Thank You for all the thing! not find "<< res.first << endl;
                            }*/
                        }
                        else{
                            filmadatype::sortPiece *curpiece = (filmadatype::sortPiece *)v.second[vi].first;
                            // 将res 插入到curpiece 的intrachain，// 将被访问的record 的数据放入内存中，即其所属的leaf，并修改其slotdata 的flag 为 true
                            gettimeofday(&lt1, NULL);

//                            curpiece->ALEX_DATA_NODE_PAYLOAD_AT(slot) = res.second;
                            curpiece->ALEX_DATA_NODE_PAYLOAD_AT(slot) = payload;
                            curpiece->ALEX_DATA_NODE_FLAG_AT(slot) = true;

                            // 将curpiece 插入到interchain 中，即mem->lru中
                            filmmem->lru->put(curpiece->ALEX_DATA_NODE_KEY_AT(0), curpiece);// 更新 interchain
                            gettimeofday(&lt2, NULL);
                            ltimeuse += (lt2.tv_sec - lt1.tv_sec) + (double) (lt2.tv_usec - lt1.tv_usec) / 1000000.0;
                            /*if (res.first != curpiece->ALEX_DATA_NODE_KEY_AT(slot)) {
                                cout << "i need You, my Lord, Thank You for all the thing! not find " << res.first << endl;
                            }*/
//                cout<< " i need You, my Lord!!! we all need You! "<< endl;
                        }
                    }
                }

                if (filmmem->check_merge(index_pageid) == true){

                    if (pageid != filmmem->inpage->pageid) {
                         // record the index_pageid, that the position in the reduced_evictedtable, which can be reused.
                        filmmem->reused_pageids.emplace_back(pageid);
                        filmmem->num_page_merge += 1;
                        filmmem->reused_index_pageid.emplace_back(index_pageid);
                        filmmem->page_merge(buf, pagedisk->recordnum, index_pageid);

//                    if (numinpage != merge_num)
//                        cout<< "Jesus, i trust in You! numinpage != merge_num" << endl;
                        filmindex->exkeynum -= filmmem->reduced_evicttable[index_pageid+1];
                        filmindex->inkeynum += filmmem->reduced_evicttable[index_pageid+1];
                        filmmem->reduced_evicttable[index_pageid+1] = 0;
                    }
//                    else{
//                         // record the index_pageid, that the position in the reduced_evictedtable, which can be reused.
//                        filmmem->reused_pageids.emplace_back(pageid);
//                        filmmem->num_page_merge += 1;
//                        filmmem->reused_index_pageid.emplace_back(index_pageid);
//                        filmmem->page_merge(filmmem->inpage->inmemdata,pagedisk->recordnum, index_pageid);
//                        filmmem->inpage->pageid +=1;
//                        filmmem->inpageid ++;
//                        filmmem->inpage->recordnum = 0;
//                        filmmem->inpage->freespace = pagedisk->pagesize;
//                        filmmem->inpage->usedspace = 0;
//
//                        filmindex->exkeynum -= filmmem->reduced_evicttable[index_pageid+1];
//                        filmindex->inkeynum += filmmem->reduced_evicttable[index_pageid+1];
//                        filmmem->reduced_evicttable[index_pageid+1] = 0;
//                    }
                }

            }
            gettimeofday(&dt2, NULL);
            dtimeuse = (dt2.tv_sec - dt1.tv_sec) + (double) (dt2.tv_usec - dt1.tv_usec) / 1000000.0;
            r_stats->disktime += dtimeuse;
            r_stats->disktime -= ltimeuse;
            r_stats->lrutime += ltimeuse;
        }
        //        cout<<"Jesus, You are my refuge! "<<endl;
        delete[] buf;
        close(fd);

        return 0;
    }


    template <class workload_t>
    int test_interleave_insert_query(unsigned int errbnd,size_t numkey,int pagesize,string dataset,double mem_threshold,double reserveMem,
                                     int recordSize,int merge_threshold, workload_t work,unsigned long actual_numkey,int datasize,double zipf,
                                     unsigned long init_num_keys,unsigned int  rebuild_threshold) {

        filmadatype* filmada = new filmadatype(0, errbnd, errbnd);
        filmadalrutype *interchain = new filmadalrutype(actual_numkey);
//        auto kkk = int (actual_numkey/600);
        filmada->leaflevel.leafpieces.reserve(int (actual_numkey/600));
        filmada->valuesize = recordSize - 1;
        filmada->buffer_ratio = work.buffer_ratio/(1-work.buffer_ratio);
        string workload = work.sample_distribution;
        std::ostringstream osse, ossr, ossd, ossm, ossp;
        osse << errbnd;
        ossp << pagesize;
        ossm << mem_threshold;
        ossr << (recordSize);
        ossd << actual_numkey;
        string diskpath = "/home/wamdm/chaohong/clionDir/updatefilm/diskpath/";
        string file_str =
                diskpath + dataset + "_totaly_rebuild_filmadapages_flattened_" + ossd.str() + "_" + ossm.str() + "_" + ossp.str() + "_" +
                ossr.str() + "_" + osse.str();
        const char *diskfile = file_str.c_str();
        int numrecord = pagesize / (recordSize);
        int bitnum = std::max(numrecord/32,1);
        filmstorage::filmmemory<key_type, key_type *, filmadatype, filmadalrutype, filmadadisk> memoryfilm(numkey,
                                                                                                           mem_threshold,
                                                                                                           filmada,
                                                                                                           interchain,bitnum);
        memoryfilm.reserveMem = reserveMem;
        memoryfilm.merge_threshold = merge_threshold;
        filmstorage::filmdisk<key_type> diskfilm(diskfile, pagesize, numrecord, recordSize, merge_threshold);
        if (numrecord <32)
            diskfilm.reduce_factor = 11;
        fstream fs;
        fs.open(diskfile, ios::in);
        if (fs) {
            fs.close();
            remove(diskfile);
        }

        int fd = -1;
        int ret = -1;
        uint64_t file_size = 2 * datasize * 1024 * 1024ULL;
        fd = open(diskfile, O_CREAT | O_RDWR, 0666);
        if (fd < 0) {
            printf("fd < 0");
            return -1;
        }

        //ret = fallocate(fd, 0, 0, file_size);
        ret = posix_fallocate(fd, 0, file_size);
        if (ret < 0) {
            printf("ret = %d, errno = %d,  %s\n", ret, errno, strerror(errno));
            return -1;
        }

        printf("fallocate create %.2fG file\n", file_size / 1024 / 1024 / 1024.0);
        close(fd);

        key_type payload[filmada->valuesize]{};
        key_type update_payload[filmada->valuesize]{1};

        //step1 bulk load the init_keys

        struct timeval oribt1, oribt2;
        double buildtimeuse;
        double initwritetime = 0.0;
        gettimeofday(&oribt1, NULL);

        memoryfilm.append(work.init_keys, payload, errbnd, errbnd);    // 执行数据bulk load ，逐个更新
        gettimeofday(&oribt2, NULL);
        buildtimeuse = (oribt2.tv_sec - oribt1.tv_sec) + (double) (oribt2.tv_usec - oribt1.tv_usec) / 1000000.0;

        // 判断bulk load 之后，是否需要transfer，并执行transfer
        pair<bool, memoryusage> transflag = memoryfilm.judgetransfer();
        unsigned int transleaves;
        double ratio;
        unsigned int leaves = transflag.second.meminfo["leaves"];
        if (memoryfilm.inpage == NULL) {
            memoryfilm.createinmempage(diskfilm.pagesize, diskfilm.recordnum);    //create a page in memory;
            memoryfilm.newpage_reduced_evicttable(0);  // memoryfilm.inpage->pageid = 0
        }

        // evitct according to lru
        while (transflag.first) {
            initwritetime += memoryfilm.filmtransfer(transflag.second.meminfo["totalusage"], &diskfilm);
            transflag = memoryfilm.judgetransfer();
        }

        ofstream savefile;
        string performance_file = "/home/wamdm/chaohong/clionDir/updatefilm/result/has_delete_filmplus_performance.txt";
        savefile.open(performance_file, ios::app);
        savefile << "method " << "FILM+flattened " << "available_memory " << mem_threshold + reserveMem << " error " << errbnd << " pagesize " << (pagesize * 8 / 1024) ;
        savefile << "k recordsize "<< recordSize << " build_time " << buildtimeuse << " dataset " << dataset << " datasize " << datasize
                 << " keynum " << numkey << " merge_threshold " << merge_threshold <<" rebuild_threshold " << rebuild_threshold << " workload " <<workload <<" ";
        savefile << "insert_ratio "<< work.insert_ratio << " ";
        savefile << "read_ratio " << work.read_ratio <<" ";
        savefile << "scan_ratio " << work.scan_ratio << " ";
        savefile << "update_ratio " << work.update_ratio<< " ";
        savefile << "scan_num " << work.scan_num<< " ";
        savefile << "buffer_ratio " << work.buffer_ratio<< " ";
        savefile << endl;
        savefile << flush;
        savefile.close();

        savefile.open(performance_file, ios::app);
        map<string, double>::iterator iter;
        savefile << "init_state " << "FILM+flattened ";
        savefile << "init_write_time " << diskfilm.initwtime / 1000000.0  << " ";
        savefile << "pages_init_num " << diskfilm.nextpageid << " ";
        savefile << "init_num_keys " << init_num_keys << " ";
        savefile << "random_seed " << work.random_seed << " ";

        for (iter = transflag.second.meminfo.begin(); iter != transflag.second.meminfo.end(); iter++)
            savefile << iter->first << " " << iter->second << " ";
        savefile << "\n";
        savefile << flush;
        savefile.close();

        // Run workload
        unsigned int i = init_num_keys;
        long long cumulative_inserts = 0, cumulative_updates = 0, retrain_updates = 0,
            cumulative_appends = 0, cumulative_ranges = 0;
        long long cumulative_lookups = 0, cumulative_deletes = 0;

        double cumulative_insert_time = 0;
        double comulative_rebuild_time = 0;
        double cumulative_lookup_time = 0, cumulative_delete_time = 0;
        double cumulative_update_time = 0, cumulative_append_time = 0;
        double cumulative_range_time = 0;
        long long cumulative_range_keys = 0;
        string lookup_distribution = workload;
        struct timeval workload_start_time, workload_run_time, insert_start_time, insert_end_time,
                lookup_start_time, lookup_end_time;
        double workload_elapsed_time = 0.0, lookup_elapsed_time = 0.0, insert_elapsed_time = 0.0,
            update_elapsed_time = 0.0, append_elapsed_time = 0.0, rebuild_elapsed_time = 0.0, delete_elapsed_time = 0.0,
            range_elapsed_time = 0;
        double time_limit = 28000;
//        double querycumulative_writetime = 0.0;
        double cumulative_query_writetime = 0.0, cumulative_delete_writetime = 0.0;
        double cumulative_insert_writetime = 0; // 记录插入process总共导致了多少 写磁盘的时间
        gettimeofday(&workload_start_time, NULL);
        unsigned int batch_no = 0;
        unsigned int total_retrain = 0;
        int total_newmodels = 0;


        access_stats query_stats;
        bool print_batch_stats = true;
        query_stats.zipffactor = zipf;
        query_stats.workload = workload;
        auto scan_num = work.scan_num;

        transflag = memoryfilm.judgetransfer(&query_stats);
        int periodV = memoryfilm.reserveMem * 1024 * 1024 / ((filmada->valuesize + 1) * sizeof(key_type) * 5);
        enum Operation {
            READ = 0, INSERT, DELETE, SCAN, UPDATE, DIFF
        };
        unsigned long int num_keys_after_batch = work.init_table_size;
        while (true) {
            batch_no++;
            if (batch_no == work.batch_num){
                cout << "i need You, my Lord, batch_no " << batch_no << endl;
                break;
            }

            // generate workload according to the work.parameters
            // generate the search key within the inserted key space

            if (i > 0) {
                double gen_duration_time = work.generate_operations(work.keys, num_keys_after_batch);

                // step1: do lookup in htap;
                double batch_compute_time = 0;
                query_stats.computetimeuse = 0.0;
                double batch_query_writetime = 0.0, bactch_lru_write_time =0, batch_delete_writetime = 0.0;
                unsigned long num_actual_lookups = work.read_ops.size();
//                unsigned long num_a = work.insert_ops.size();
                cumulative_range_keys += num_actual_lookups;
                if (num_actual_lookups) {
                    int qi = 0;
                    gettimeofday(&lookup_start_time, NULL);
                    for (auto &op_key: work.read_ops) {
//                        if ( op_key == 17505691502760){
//                            cout << "Jesus, where this key go?" << endl;
//                        }
                        pair<key_type, key_type *> res = memoryfilm.get_key(op_key, &query_stats, &diskfilm);
                        if (query_stats.disknum > 0 && qi++ % periodV == 0) {
//                            gettimeofday(&insert_start_time, NULL);
                            transflag = memoryfilm.judgetransfer(&query_stats);
                            while (transflag.first) {
                                batch_query_writetime += runtimeevictkeytopage2(&memoryfilm, transflag.second.meminfo["totalusage"],
                                                                                filmada, &diskfilm, &query_stats);
                                transflag = memoryfilm.judgetransfer(&query_stats);
                            }
//                            gettimeofday(&insert_end_time, NULL);
//                            bactch_lru_write_time += (insert_end_time.tv_sec - insert_start_time.tv_sec) +
//                                                  (double) (insert_end_time.tv_usec - insert_start_time.tv_usec) / 1000000.0;

                        }
//                        cout << "this is read, thank You, my Lord" << endl;
                    }

                    // record the start time of lookup
                    gettimeofday(&lookup_end_time, NULL);
                    lookup_elapsed_time = (lookup_end_time.tv_sec - lookup_start_time.tv_sec) +
                                          (double) (lookup_end_time.tv_usec - lookup_start_time.tv_usec) / 1000000.0;
                    lookup_elapsed_time -= query_stats.computetimeuse;
                    lookup_elapsed_time += batch_query_writetime;
                    cumulative_lookup_time += lookup_elapsed_time;
                    cumulative_lookups += num_actual_lookups;
                    cumulative_query_writetime += batch_query_writetime;
                    vector<key_type>().swap(work.read_ops);
                    batch_compute_time += query_stats.computetimeuse;
                    query_stats.computetimeuse = 0;
                    for(iter = transflag.second.meminfo.begin(); iter !=  transflag.second.meminfo.end(); iter++)
                        cout<<iter->first<<" "<<iter->second<<" *** ";
                    cout<<endl;
                }

//              do deletes
                unsigned long num_actual_deletes = work.delete_ops.size();
                cumulative_range_keys += num_actual_deletes;
                if (num_actual_deletes){
                    int qi = 0;
                    gettimeofday(&lookup_start_time, NULL);
                    for (auto &op_key: work.delete_ops) {
//                        if ( op_key ==  16147608666968){
//                            cout << "Jesus, where this key go?" << endl;
//                        }
//                        if (qi++ < 15)
//                            cout << op_key << endl;
                        bool del_res = memoryfilm.delete_key(op_key, &diskfilm);

                        // 判断是否需要在delete 的过程中 执行transfer, 也就是eviction process: 在delete 的过程中，不需要
//                        if (query_stats.disknum > 0 && qi++ % periodV == 0) {
//                            transflag = memoryfilm.judgetransfer(&query_stats);
//                            while (transflag.first) {
//                                batch_delete_writetime += runtimeevictkeytopage2(&memoryfilm, transflag.second.meminfo["totalusage"],
//                                                                                filmada, &diskfilm, &query_stats);
//                                transflag = memoryfilm.judgetransfer(&query_stats);
//                            }
//                        }
//                        cout << "this is delete, thank You, my Lord" << endl;
                    }

                    // record the start time of lookup
                    gettimeofday(&lookup_end_time, NULL);
                    delete_elapsed_time = (lookup_end_time.tv_sec - lookup_start_time.tv_sec) +
                                          (double) (lookup_end_time.tv_usec - lookup_start_time.tv_usec) / 1000000.0;
                    delete_elapsed_time -= query_stats.computetimeuse;
                    delete_elapsed_time += batch_delete_writetime;
                    cumulative_delete_time += delete_elapsed_time;
                    cumulative_deletes += num_actual_deletes;
                    cumulative_delete_writetime += batch_delete_writetime;
                    vector<key_type>().swap(work.delete_ops);
                    batch_compute_time += query_stats.computetimeuse;
                    query_stats.computetimeuse = 0;
                    for(iter = transflag.second.meminfo.begin(); iter !=  transflag.second.meminfo.end(); iter++)
                        cout<<iter->first<<" "<<iter->second<<" *** ";
                    cout<<endl;
                }

                    // step2: Do inserts, with htap
                unsigned last_retrain = filmada->retrain;
                int last_newmodels = filmada->newmodels;
                double batch_insert_writetime = 0.0;   // sec
                unsigned long int num_actual_inserts = work.insert_ops.size();
                int batch_newmodels = 0;
                unsigned batch_retrain = 0;
                cumulative_range_keys += num_actual_inserts;
                if (num_actual_inserts){
//                    gettimeofday(&insert_start_time, NULL);   //
                    auto start_time = std::chrono::high_resolution_clock::now();
//                   filmadatypeinsert_random(work.insert_ops, payload, interchain);
                    for (auto &op_key: work.insert_ops){
//                        if (op_key == 12368769653296 ) { //  15945001399101
//                             cout << "my Lord, i need You! thank You, for ever and ever ~~" << endl;
//                        }
//                       filmada->insert_one(op_key,payload,interchain);
                        filmada->insert(op_key,payload,interchain);
                    }
                    filmada->inkeynum += num_actual_inserts;
//                    cout << "i need You, my Lord" << endl;
                    // used in insert_one procedure, update the last piece
                    /*
                   filmada->root =filmada->innerlevels.back()->innerpieces[0];
                   filmada->inkeynum += num_actual_inserts;
                    auto a =filmada->leaflevel.opt->get_segment(filmada->m_tailleaf->endkey);
                    if (filmada->leaflevel.opt->points_in_hull< 2){
                        a.first =filmada->m_tailleaf->endkey;
                        a.last =filmada->m_tailleaf->endkey;
                    }
                   filmada->m_tailleaf->update(a);
                    */
                    insert_elapsed_time =
                    std::chrono::duration_cast<std::chrono::nanoseconds>(
                            std::chrono::high_resolution_clock::now() - start_time)
                            .count();

                    // step3: judge whether to evict data, after doing inserts
                    transflag = memoryfilm.judgetransfer(&query_stats);
                    query_stats.computetimeuse = 0.0;
                    while (transflag.first) {
                        batch_insert_writetime += runtimeevictkeytopage2(&memoryfilm, transflag.second.meminfo["totalusage"],
                                                                         filmada, &diskfilm, &query_stats);
                        transflag = memoryfilm.judgetransfer(&query_stats);
                    }

                    vector<key_type>().swap(work.insert_ops);

                    batch_newmodels =filmada->newmodels - last_newmodels;
                    batch_retrain =filmada->retrain - last_retrain;
                    total_retrain += batch_retrain;
                    total_newmodels += batch_newmodels;

                    if (filmada->retrain >= rebuild_threshold){

                        auto start_time = std::chrono::high_resolution_clock::now();
                        // rebuild internal levels
                        filmada->internal_level_rebuild(errbnd,interchain);

                        rebuild_elapsed_time = std::chrono::duration_cast<std::chrono::nanoseconds>(
                                std::chrono::high_resolution_clock::now() - start_time)
                                .count();
//                        filmada->newmodels = 0;
                        filmada->retrain = 0;
                        filmada->rebuild_inner += 1;
                        filmada->key_in_buff = 0;
                        filmada->reserve_buff = 0;
                        filmada->fill_factor = 0.0;
                    }

                    // 将写出磁盘的时间算到插入时间里，更新插入的时间
                    insert_elapsed_time /= 1e9;
                    rebuild_elapsed_time /= 1e9;
                    insert_elapsed_time += rebuild_elapsed_time;
                    comulative_rebuild_time += rebuild_elapsed_time;

                    insert_elapsed_time += batch_insert_writetime;
//                insert_elapsed_time -= query_stats.computetimeuse;

                    cumulative_inserts += num_actual_inserts;
                    num_keys_after_batch += num_actual_inserts;

                    cumulative_insert_time += insert_elapsed_time;
                    cumulative_insert_writetime += batch_insert_writetime;
                }


                // step 4 do update, change the payload
                unsigned long int num_actual_updates = work.update_ops.size();
                cumulative_range_keys += num_actual_updates;
                if (num_actual_updates){
                    gettimeofday(&insert_start_time, NULL);   //
                    for (auto &op_key:work.update_ops) {
//                        if (op_key == 291660467280 ) {
//                            cout << "my Lord, i need You! thank You, for ever and ever ~~" << endl;
//                        }
                        // 先找到key，再更新
                        pair<key_type, key_type *> res = memoryfilm.update_key(op_key, update_payload, &query_stats,diskfilm);
//                        cout << "this is update, thank You, my Lord" << endl;
                    }
                    gettimeofday(&insert_end_time, NULL);
                    update_elapsed_time = (insert_end_time.tv_sec - insert_start_time.tv_sec) +
                                          (double) (insert_end_time.tv_usec - insert_start_time.tv_usec) / 1000000.0;
                    cumulative_update_time += update_elapsed_time;
                    cumulative_updates += num_actual_updates;
//                cumulative_appends += alen;
                    cumulative_insert_writetime += batch_insert_writetime;

                    vector<key_type>().swap(work.update_ops);
                }
//                step 5 do range scan
                batch_query_writetime = 0;
                bactch_lru_write_time = 0;
                unsigned long int num_actual_ranges = work.scan_ops.size();
                query_stats.computetimeuse = 0.0;
                if (num_actual_ranges){
                    gettimeofday(&insert_start_time, NULL);   //
//                    for (auto &op_key:work.scan_ops) {
////                        pair<key_type, key_type *> res = memoryfilm.get_scan(op_key, &query_stats,diskfilm);
//                        cout << "this is scan, thank You, my Lord" << endl;
//                    }
                    for (unsigned int ri = 0; ri < num_actual_ranges; ri++) {
//                        if (batch_no == 8 && ri == 32768){
//                            COUT_THIS("Jesus, i believe, please help my unbelief.");
//                        }
                         vector<pair<key_type, key_type *>> rres;
                        key_type firstkey = work.range_ops[ri][0];
                        key_type lastkey = work.range_ops[ri][1];
                        auto index_res = memoryfilm.get_range(firstkey,lastkey, &query_stats,diskfilm);

//                        auto index_res = memoryfilm.get_scan(firstkey,scan_num, &query_stats,diskfilm);
                        pre_dict *dict = new pre_dict();
                        htap_range_prepass(&memoryfilm, &index_res, &rres, dict, filmada, &query_stats,diskfilm.reduce_factor);

                        // determine whether to access disk
                        if (dict->pagenum == 0) {  // the request data  数据都在内存
                            //          cout<< "Jesus, happy birthday!" << endl;
                            query_stats.memnum += 1;
                            // 遍历rres，释放数组
                        } else if (rres.empty()) { //数据都在磁盘
                            range_read_from_disk(&rres, dict, &memoryfilm, filmada, &diskfilm,
                                                 &query_stats,payload);// read data from disk according to 根据prepass 的信息，从磁盘中读数据

                            query_stats.disknum += 1;
                            query_stats.diskpagenum += dict->pagenum;
                            //cout << "i want in Your heart, my Lord!" << endl;
                        } else {   //the request data are in both memory and disk
                            query_stats.crossnum += 1;
                            query_stats.crosspagenum += dict->pagenum;
                            //                cout<< "Jesus, sister nana needs You!"<<endl;
                            range_read_from_disk(&rres, dict, &memoryfilm, filmada, &diskfilm, &query_stats,payload);
                        }
                        delete dict;
                        cumulative_range_keys += rres.size();
                        vector<pair<key_type, key_type *>>().swap(rres) ;
//                        if (rres[0].first < firstkey)
//                            COUT_THIS("Jesus, please come");

                        if (query_stats.disknum > 0 && ri % periodV == 0) {
//                            gettimeofday(&lookup_start_time, NULL);
                            transflag = memoryfilm.judgetransfer(&query_stats);
                            while (transflag.first) {
                                batch_query_writetime += runtimeevictkeytopage2(&memoryfilm, transflag.second.meminfo["totalusage"],
                                                                                filmada, &diskfilm, &query_stats);
                                transflag = memoryfilm.judgetransfer(&query_stats);
                            }
//                            gettimeofday(&lookup_end_time, NULL);
//                            bactch_lru_write_time += (lookup_end_time.tv_sec - lookup_start_time.tv_sec) +
//                                                     (double) (lookup_end_time.tv_usec - lookup_start_time.tv_usec) / 1000000.0;
                        }
                    }
                    gettimeofday(&insert_end_time, NULL);
                    range_elapsed_time = (insert_end_time.tv_sec - insert_start_time.tv_sec) +
                                          (double) (insert_end_time.tv_usec - insert_start_time.tv_usec) / 1000000.0;

                    range_elapsed_time -= query_stats.computetimeuse;
                    range_elapsed_time += batch_query_writetime;
                    cumulative_query_writetime += batch_query_writetime;
                    cumulative_range_time += range_elapsed_time;
                    cumulative_ranges += num_actual_ranges;
                   for (unsigned int ri = 0; ri < num_actual_ranges; ri++)
                        delete[] work.range_ops[ri];
                    delete[] work.range_ops;
                    work.range_ops = nullptr;
                    vector<key_type>().swap(work.scan_ops);
                    batch_compute_time += query_stats.computetimeuse;
                    query_stats.computetimeuse = 0;
                }



                // the total time
                gettimeofday(&workload_run_time, NULL);
                workload_elapsed_time = (workload_run_time.tv_sec - workload_start_time.tv_sec) +
                                        (double) (workload_run_time.tv_usec - workload_start_time.tv_usec) / 1000000.0;
                workload_elapsed_time -= gen_duration_time;
                workload_elapsed_time -= batch_compute_time ;


                struct timeval fill_start_time,fill_end_time;
//                gettimeofday(&fill_start_time, NULL);
////                if (batch_no % 10 == 0)
////                    filmada->compute_fill_factor();
//                gettimeofday(&fill_end_time, NULL);
//                double fill_elapsed_time = (fill_end_time.tv_sec - fill_start_time.tv_sec) +
//                                               (double) (fill_end_time.tv_usec - fill_start_time.tv_usec) / 1000000.0;
//                workload_elapsed_time -= fill_elapsed_time;


                // step 6: print 输出  batch_no ==1 || batch_no%40 == 0
                if (true) {
                    savefile.open(performance_file, ios::app);
                    map<string, double>::iterator iter;
//                savefile << "_ ";
                    savefile << "disk_write_time " << diskfilm.initwtime / 1000000.0  << " ";
                    savefile << "used_page_num " << (diskfilm.nextpageid-memoryfilm.reused_pageids.size()) << " ";
                    savefile << "method " << "FILM+flattened ";
                    for (iter = transflag.second.meminfo.begin(); iter != transflag.second.meminfo.end(); iter++)
                        savefile << iter->first << " " << iter->second << " ";
//                savefile << "\n";
                    savefile << flush;
                    savefile.close();
                    query_stats.print_stats();
                    int num_batch_operations = num_actual_lookups + num_actual_inserts +
                                        num_actual_updates + num_actual_ranges + num_actual_deletes;
                    double batch_time = lookup_elapsed_time + insert_elapsed_time + delete_elapsed_time +
                                        update_elapsed_time + range_elapsed_time;
                    long long cumulative_operations = cumulative_lookups + cumulative_inserts + cumulative_deletes +
                                    cumulative_updates + cumulative_ranges;
                    double cumulative_time = cumulative_lookup_time + cumulative_insert_time + cumulative_delete_time +
                            cumulative_update_time + cumulative_range_time;

                    std::cout << "Batch " << batch_no
                              << ", cumulative_ops: " << cumulative_operations
                              << "\n\tbatch_throughput:\t"
                              << num_actual_lookups / lookup_elapsed_time
                              << " lookups/sec,\t"
                              << num_actual_inserts / insert_elapsed_time
                              << " inserts/sec,\t"
                                << num_actual_deletes / delete_elapsed_time
                                << " deletes/sec,\t"
                              << (num_actual_updates) / (update_elapsed_time)
                              << " updates/sec,\t"
                            << (num_actual_ranges) / (range_elapsed_time)
                            << " ranges/sec,\t"
                            << num_batch_operations / batch_time
                              << " ops/sec "
                              << "\n\tcumulative_throughput:\t"
                              << cumulative_lookups / cumulative_lookup_time
                              << " lookups/sec,\t"
                              << cumulative_inserts / cumulative_insert_time
                              << " inserts/sec,\t"
                            << cumulative_deletes / cumulative_delete_time
                            << " deletes/sec,\t"
                              << cumulative_updates / cumulative_update_time
                              << " updates/sec,\t"
                            << cumulative_ranges / cumulative_range_time
                            << " ranges/sec,\t"
                              << cumulative_operations / cumulative_time << " ops/sec "
                            << cumulative_range_keys / cumulative_time << " cum_keys/sec "
                              << cumulative_query_writetime << " cumulative_query_writetime "
                              << cumulative_insert_writetime << " cumulative_insert_writetime "
                              << rebuild_elapsed_time << " rebuild_elapsed_time "
                                << delete_elapsed_time << " delete_elapsed_time ";

                    std::cout << "\n\tcumulative_execution_1:\t"
                              << ((filmada->inkeynum +filmada->exkeynum)/double(1024 * 1024)) * (filmada->valuesize + 1) * 8
                              << " G datasize, \t"
                              <<filmada->inkeynum +filmada->exkeynum << " #cumulative_keys,\t"
                              <<filmada->inkeynum << " #inmem_keys,\t"
                              <<filmada->exkeynum << " #exmem_keys"
                              << "\n\tcumulative_execution_1:\t"
                              << cumulative_updates << " #cumulative_updates,\t" << cumulative_update_time
                              << " cumulative_update_time,\t"
                              << cumulative_inserts << " #cumulative_inserts,\t" << cumulative_insert_time
                              << " cumulative_insert_time,\t"
                            << cumulative_deletes << " #cumulative_deletes,\t" << cumulative_delete_time
                            << " cumulative_delete_time,\t"
                              << cumulative_lookups << " #cumulative_lookups,\t" << cumulative_lookup_time
                              << " cumulative_lookup_time,\t"
                            << cumulative_ranges << " #cumulative_ranges,\t" << cumulative_range_time
                            << " cumulative_range_time,\t"
                             << total_retrain << " cumulative_retrain,\t "
                              <<  comulative_rebuild_time << " comulative_rebuild_time,\t"
                              << (diskfilm.initwtime / 1000000.0) + query_stats.wdisktime << " cumulative_writetime "
                              << cumulative_insert_writetime << " cumulative_insert_writetime "
                              << cumulative_delete_time << " cumulative_delete_time ";

                    std::cout << "\n\tfilm_update_1:\t"
                              << rebuild_elapsed_time << " rebuild_elapsed_time "
                              << memoryfilm.num_page_merge << " num_page_merge "
                              << memoryfilm.merge_timeuse << " merge_timeuse "
                              << total_retrain << " times_of_leaves_retrain "
                              << total_newmodels << " cumulative_new_models "
                              << batch_retrain<< " batch_times_of_leaves_retrain "
                              << batch_newmodels << " batch_updated_new_models "
                            <<filmada->reserve_buff << " reserve_buff "
                            <<filmada->key_in_buff << " key_in_buff "
                            <<filmada->fill_factor << " fill_factor ";

                    std::cout << std::endl;

                    savefile.open(performance_file,ios::app);

                    savefile << batch_no << " Batch " <<
                             cumulative_operations << " cumulative_ops "
                             << "_" << " batch_throughput "
                             << num_actual_lookups / lookup_elapsed_time << " lookups/sec "
                             << num_actual_inserts / insert_elapsed_time << " inserts/sec "
                            << num_actual_deletes / delete_elapsed_time << " deletes/sec "
                             << num_actual_updates / update_elapsed_time << " updates/sec "
                            << num_actual_ranges / range_elapsed_time << " ranges/sec "
                             << num_batch_operations / batch_time << " ops/sec "
                             << batch_query_writetime << " batch_query_writetime "
                             << batch_insert_writetime << " batch_insert_writetime "
                             << "_" << " cumulative_throughput "
                             << cumulative_lookups / cumulative_lookup_time << " cum_lookups/sec "
                             << cumulative_inserts / cumulative_insert_time << " cum_inserts/sec "
                            << cumulative_deletes / cumulative_delete_time << " cum_deletes/sec "
                             << cumulative_updates / cumulative_update_time << " cum_updates/sec "
                            << cumulative_ranges / cumulative_range_time << " cum_ranges/sec "
                             << cumulative_operations / cumulative_time << " cum_ops/sec "
                            << cumulative_range_keys / cumulative_time << " cum_keys/sec "
                             << cumulative_query_writetime << " cumulative_query_writetime "
                             << cumulative_insert_writetime << " cumulative_insert_writetime "
                             << rebuild_elapsed_time << " rebuild_elapsed_time "
                             <<filmada->rebuild_inner << " num_of_inner_rebuild "
                            << total_retrain << " cumulative_retrain "
                            << total_newmodels << " cumulative_new_models "
                            << total_retrain << " cumulative_retrain "
                            <<  comulative_rebuild_time << " comulative_rebuild_time"
                            << batch_retrain<< " batch_times_of_leaves_retrain "
                            << batch_newmodels << " batch_updated_new_models "
                            <<filmada->maxbsize  << " maxbuffersize "
                            <<filmada->minbsize  << " minbuffersize "
                            <<filmada->reserve_buff << " reserve_buff "
                            <<filmada->key_in_buff << " key_in_buff "
                            <<filmada->fill_factor << " buff_fill_factor ";

                    auto allkeynum =filmada->inkeynum +filmada->exkeynum;
                    savefile << (filmada->inkeynum +filmada->exkeynum) * (filmada->valuesize + 1) * 8 / double(1024 * 1024)
                             << "G datasize "
                             << allkeynum << " #cumulative_keys "
                             << allkeynum/double(allkeynum+filmada->reserve_buff) << " global_fill_factor "
                             <<filmada->inkeynum << " #inmem_keys "
                             <<filmada->exkeynum << " #exmem_keys "
                             << cumulative_lookups << " #cumulative_lookups " <<
                             cumulative_lookup_time  << " cumulative_lookup_time "
                            << cumulative_inserts << " #cumulative_inserts "
                            << cumulative_insert_time << " cumulative_insert_time "
                            << cumulative_deletes << " #cumulative_deletes "
                            << cumulative_delete_time << " cumulative_delete_time "
                            << cumulative_updates << " #cumulative_updates "
                            << cumulative_update_time<< " cumulative_update_time "
                            << cumulative_ranges << " #cumulative_ranges "
                            << cumulative_range_keys  << " cumulative_range_keys "
                            << cumulative_range_time<< " cumulative_range_time "
                             << (diskfilm.initwtime / 1000000.0) + query_stats.wdisktime << " cumulative_writetime "
                             << memoryfilm.num_page_merge << " num_page_merge "
                             << memoryfilm.merge_timeuse << " merge_timeuse ";

                    savefile << "\n";
                    savefile << flush;
                    savefile.close();

                }


                // Check for workload end conditions
                if (workload_elapsed_time > time_limit || (filmada->inkeynum +filmada->exkeynum) >= actual_numkey) {
                    break;
                }

//                if (num_actual_inserts < num_inserts_per_batch) {
//                    // End if we have inserted all keys in a workload with inserts
//                    break;
//                }
            }
            cout << "my lovely Lord, finished the interleave inserts and queries of batch " << batch_no << endl;

        }
        for (iter = transflag.second.meminfo.begin(); iter != transflag.second.meminfo.end(); iter++)
            savefile << iter->first << " " << iter->second << " ";
        savefile << "\n";
        savefile << flush;
        savefile.close();
        vector<unsigned int>().swap(memoryfilm.reduced_evicttable);
        fs.open(diskfile,ios::in);
        delete filmada;
        if (fs){
            fs.close();
            remove(diskfile);
        }
        malloc_trim(0);
        return 0;
    }



}

#endif //PLUSFILM_FILM_PLUS_H

/*
 * search_one_leaf  consider the regular data update
predicted_pos = (*a)(key);

int slotlo = PGM_SUB_EPS2(predicted_pos, ErrorRecursive + 1,(a->slotkey.size()-1)) ;  //,(--itlevel)->size()
int slothi = PGM_ADD_EPS(predicted_pos, ErrorRecursive,(a->slotkey.size()-1));
int resslot = slotlo;
for (; resslot<slothi;resslot++){
if (key <= a->slotkey[resslot])
{
break;
}
}
result_find index_res;
if (key != a->slotkey[resslot]){
index_res= result_find(false,a->locbitmap[resslot],a,resslot);   // not find, then, insert into bufferpiece,
index_res.leafslot = search_pos;
}
else{
index_res = result_find(true,a->locbitmap[resslot],a,resslot);   // regular data
index_res.leafslot = search_pos;
}

return index_res;
 */

/*
forceinline result_find search_one(const key_type key) {
    // 从 root 开始，即，从root开始的结尾开始

    auto itlevel = innerlevels.end() -1; // root level 的 iterator
    // internal level use exponential search, the parameters used in exponential_search are: predict_pos, search key
//
//     * 1. get the predicted pos
//     * 2. get the [l,r] used in binary search or, get the [l,r]
//     * 3. use the binary search in [l,r] to get the final position
//
    auto predicted_pos = (*root)(key);
    auto childpieces = (*(--itlevel))->innerpieces;
//                auto search_pos = std::max(exponential_search_upper_bound(childpieces,predicted_pos, key) -1,0);
    auto search_pos =exponential_search_upper_bound(childpieces,predicted_pos, key) -1;

    // 判断 找到的inner node 是否满足条件
    auto curit = (*(itlevel))->innerpieces.begin() + search_pos;
    Innerpiece* b = *curit;

    for (; --itlevel >= innerlevels.begin();){
        predicted_pos = (*b)(key);
        childpieces = (*(itlevel))->innerpieces;
        search_pos = exponential_search_upper_bound(childpieces,predicted_pos, key) -1;

        // 判断 找到的inner node 是否满足条件
        curit = (*(itlevel))->innerpieces.begin() + search_pos;
    }
    // find the leaf level

    b = *curit;
    predicted_pos = (*b)(key);

    search_pos = exponential_search_upper_bound(leaflevel.leafpieces,predicted_pos, key) -1;

    // 预测位置和最终位置的距离
    auto curleaf = leaflevel.leafpieces.begin() + search_pos;

    Leafpiece* a = *curleaf;
    predicted_pos = (*a)(key);
    // 确定 该leaf 的哪一个链 leaf piece
    while ( a->buffer_piece->buffer_piece && key > a->endkey) {
        a = a->buffer_piece;
    }
    if (a->buffer_piece == NULL) {
//                    gettimeofday(&ct1, NULL);
        sortPiece *buffer = (sortPiece *) a;
        auto resslot = predicted_pos * this->buffer_ratio;
        auto model_res = buffer->find_key(key, resslot);

//                if (buffer->ALEX_DATA_NODE_KEY_AT(model_res) == key) {
        result_find index_res = result_find(false, buffer->ALEX_DATA_NODE_FLAG_AT(model_res), buffer,
                                            model_res);   // sort_list data
        return index_res;


    }
    else {
        predicted_pos = (*a)(key);

        int slotlo = PGM_SUB_EPS(predicted_pos, ErrorRecursive + 1);
        slotlo = std::min<int>(slotlo,(a->slotkey.size() - 1));
        int slothi = PGM_ADD_EPS(predicted_pos, ErrorRecursive, (a->slotkey.size() - 1));
        int resslot = slotlo;
        for (; resslot < slothi; resslot++) {
            if (key <= a->slotkey[resslot]) {
                break;
            }
        }
//                    gettimeofday(&ct1, NULL);
        if (key != a->slotkey[resslot]) {
            // find in sort_list
            while (a->buffer_piece->buffer_piece != NULL) {
                a = a->buffer_piece;
            }
            predicted_pos = (*(Leafpiece*)(*curleaf))(key);
            sortPiece *buffer = (sortPiece *) a->buffer_piece;
            resslot = predicted_pos * this->buffer_ratio;
            auto model_res = buffer->find_key(key, resslot);
//                    if (buffer->ALEX_DATA_NODE_KEY_AT(model_res) == key) {
            result_find index_res = result_find(false, buffer->ALEX_DATA_NODE_FLAG_AT(model_res), buffer,
                                                model_res);   // sort_list data
            return index_res;

        }

        if ( a->delbitmap[resslot]){
            result_find index_res = result_find(true, a->locbitmap[resslot], a, -1);   // regular data
            return index_res;
        }
        result_find index_res = result_find(true, a->locbitmap[resslot], a, resslot);   // regular data
        return index_res;
    }

}
*/