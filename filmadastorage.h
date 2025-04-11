

//
// Created by chaochao on 2021/12/23.
//

#ifndef EXPERIMENTCC12_FILMADASTORAGE_H
#define EXPERIMENTCC12_FILMADASTORAGE_H

#include <fcntl.h>
#include <unistd.h>
#include <memory.h>
//#define _GNU_SOURCE
#include <vector>
#include "filmadalru.h"
#include "film_plus.h"


using namespace std;
typedef unsigned lruOff_type;
typedef unsigned short int pageOff_type;
typedef unsigned int pageid_type;
namespace filmstorage {
    typedef std::map<std::string,double> infomap;
    template <class key_type,class data_type,class index_type, class lru_type,class disk_type>
    class filmmemory{
    public:
        int totalnum;
        pageid_type inpageid = 0;
        double threshold;
        double reserveMem;
        int merge_threshold;
        double merge_timeuse = 0.0;
        unsigned int num_page_merge = 0;
        index_type *index;  //
        lru_type *lru;
        vector<pair<pageid_type ,pageOff_type>> evicttable;  // pageid , offset in page
        vector<pageid_type> reduced_evicttable;
        vector<pageid_type> reused_pageids;
        vector<unsigned long> reused_index_pageid;
        // reduced_evicttable: pos_0 page_id;pos_1 num_in_page; pos2-257 bitmap; pos258-4353 offsets
        // a page_information_block has 4353个unsigned int, 共 17416 个字节，
        // 按照unsigned short 来计数，第一个block的index范围为0-8705，每8706个index代表一个block
        // 每一个offset 占用unsigned int 类型的高16位或低16位。
        // the information of record on disk, is the index (indexed as unsigned short) in the reduced_evictedtable,
        // thus, the page_id is [0]~(0,1), the num in page is in reduced_evict[1]~(2,3),
        // the bitmap is reuced_evicttable[2-257]
        // the offset is in one of the poss (258-4353)
        // the information of the second page, the offset is in one of the poss (8706+4-17411)
        // the information of the i-th page (i=[0,page_num]), the offsets are in the poss [8706*i+4,8706*i-1]
        // 根据offset 要读取page 和该page 中的record 时，for example, the index is i_pos = 8709. i_pos/8706 得到整数位和余数位分别为 1 和 4
        // 读取reduceevicttable[4353*1+4/2], pageid 为4353*1;
        // track_index = 8706*i + k, (k = 2x+c, c = 0 or 1) the page id is in reduceevictable[4353*i], the offset is in reduceevicttable[4353*i+(k/2)+c]
        vector<key_type> evictkey;   // the evicted key, used in debug
        vector<pair<pageid_type ,pageOff_type>*> evictPoss;
        vector<unsigned int> bitmap;
        int bitsnum;


        filmmemory( unsigned long int numkey,double Threshold, index_type *filmtree , lru_type *LRU, int bitnum){
            totalnum = 0;
            threshold = Threshold;
//            evicttable.reserve(numkey);
            reduced_evicttable.reserve(numkey);
            index = filmtree;
            lru = LRU;
            for (int i = 0; i < bitnum;i ++){
                bitmap.emplace_back(0);
            }
            bitsnum = bitnum;
            bitmap.resize(bitnum);
        }

        void insert(const std::vector<key_type> keys,key_type* payload,const int error,const int error_recursize ){
            index->build(keys,payload,error,error_recursize,lru);
        }

        void append(const std::vector<key_type> keys,key_type* payload,const int error,const int error_recursize ){
            index->update_append(keys,payload,error,error_recursize,lru);
        }

        void update(const std::vector<key_type> keys,key_type* payload){
            index->update_random(keys,payload,lru);    // out-of-order insertion
        }

        std::tuple<pageid_type,pageOff_type, unsigned long  > query_flat_structure(unsigned long short_reduced_evict, double reduce_factor){
            auto int_reduced_evict = short_reduced_evict/2.0;
            unsigned long index_pageid = floor(int_reduced_evict/reduce_factor)*reduce_factor;
            pageid_type pageid = reduced_evicttable[index_pageid];
            pageOff_type offset = 0;
            auto fraction = modf(int_reduced_evict,&int_reduced_evict);
            auto int_high_low = reduced_evicttable[int_reduced_evict];
            // 根据reduced_evict 算出（读出） page id，算出读出offset
            if (fraction){
                offset = int_high_low&0xFFFF;
            }
            else
                offset =  int_high_low>>16 ;   // 高16位
            return {pageid, offset,index_pageid};
        }


        template< typename stat_type>
        pair<key_type, key_type *> get_key(key_type key, stat_type query_stats, disk_type* diskfilm){
//            struct timeval tqt1, tqt2, xt1, xt2, plt1, plt2, prdt1, prdt2;
//            double pltimeuse, prdtimeuse, pwdtimeuse, timeuse, xtimeuse;  // 读磁盘的时间，写磁盘的时间
            pair<key_type, key_type *> res;
//            gettimeofday(&xt1, NULL);
//            auto index_res1 = index->search_one(key, query_stats);
            auto index_res = index->search_one(key);
//            gettimeofday(&xt2, NULL);
//            xtimeuse = (xt2.tv_sec - xt1.tv_sec) + (double) (xt2.tv_usec - xt1.tv_usec) / 1000000.0;
//            query_stats->xlookuptime += xtimeuse;
            if (index_res.slot == -1){
                return res;
            }
            else{
                if (index_res.find == false) {   // 从 sort_piece 中 读取数据
                    typename index_type::sortPiece* findpiece = (typename index_type::sortPiece*) (index_res.findleaf);
//                    auto res_key = findpiece->ALEX_DATA_NODE_KEY_AT(index_res.slot);
//                    if (res_key != key)   // 说明已经删除了key
//                    {
//                        COUT_THIS("res_key != key~ Lord, please come, You are my refuge from ever to ever, ")
//                        return res;
//                    }

                    if (index_res.flags) {
                        query_stats->memnum += 1;
                        auto finddata = (key_type *) (findpiece->ALEX_DATA_NODE_PAYLOAD_AT(index_res.slot));
                        res.first = key;
                        res.second = finddata;
                        // 更新 interchain
                        lru->put(findpiece->ALEX_DATA_NODE_KEY_AT(0), findpiece);
                    }
                    else {
                        query_stats->disknum += 1;
                        query_stats->diskpagenum += 1;

                        typename index_type::sortPiece* findpiece = (typename index_type::sortPiece*) (index_res.findleaf);

//                        auto res_key = findpiece->ALEX_DATA_NODE_KEY_AT(index_res.slot);
//                        if (res_key != key)
//                            COUT_THIS("res_key != key~ Lord, please come, You are my refuge from ever to ever, ")

                        auto short_reduced_evict = (unsigned long) findpiece->ALEX_DATA_NODE_PAYLOAD_AT(index_res.slot);  //原本这里可以直接是 writeevict 指向的是 pageid and offset
                        // 在FILM+ flattened structures 中， 为进一步降低 memory usage，改成了flattened structure 的结构
                        auto [pageid, offset,index_pageid] = query_flat_structure(short_reduced_evict,diskfilm->reduce_factor);

                        // 将 querykey merge到memory，修改 bitmap，修改num in page
                        reduced_evicttable[index_pageid+1] -= 1;

                        if (pageid == inpage->pageid) {
                            res.second = new key_type[index->valuesize];
                            key_type *inmemdata = inpage->inmemdata;
                            res.first = inmemdata[offset];
                            for (int ki = 0; ki < index->valuesize; ki++) {
                                res.second[ki] = inmemdata[offset + 1 + ki];
                            }
                            auto recordnum = offset/diskfilm->recordsize;
                            this->unset_bit(recordnum, index_pageid);
                        }
                        else {
                            res = diskfilm->odirectreadfromdisk(pageid,offset,this,index_pageid);    // if readfromdisk indicating doesn't use o_direct;
                        }

                        findpiece->ALEX_DATA_NODE_PAYLOAD_AT(index_res.slot) = res.second;
                        findpiece->ALEX_DATA_NODE_FLAG_AT(index_res.slot) = true;
                        lru->put(findpiece->ALEX_DATA_NODE_KEY_AT(0), findpiece);// 更新 interchain

                        index->inkeynum++;
                        index->exkeynum--;
                    }
                }
                else {
                    if (index_res.flags) {// 从内存中读数据
                        query_stats->memnum += 1;
                        auto finddata = (adalru::Node<lruOff_type, key_type *> *) index_res.findleaf->slotdata[index_res.slot];
                        res.first = key;
                        bool del_flag = index_res.findleaf->delbitmap[index_res.slot];
//                        if (del_flag)
//                            return res;
//                        else{
                        res.second = finddata->value;
                        index_res.findleaf->intrachain.moveTohead(finddata);// update intrachain
                        // 更新 interchain
                        lru->put(index_res.findleaf->startkey, index_res.findleaf);
//                        }
                    }
                    else { //从磁盘读数据
                        query_stats->disknum += 1;
                        query_stats->diskpagenum += 1;

                        auto short_reduced_evict = (unsigned long) index_res.findleaf->slotdata[index_res.slot];   //原本这一步可以直接 指向的是 pageid and offset

                        auto [pageid, pageoff,index_pageid] = query_flat_structure(short_reduced_evict,diskfilm->reduce_factor);
//                    gettimeofday(&xt1, NULL);

                        // 将 querykey merge到memory：修改num in page， 修改 records bitmap (在读磁盘的时候进行了修改 memoryfilm->unset_bit())
                        this->reduced_evicttable[index_pageid+1] -= 1;

//                    gettimeofday(&prdt1, NULL);
                        if (pageid == inpage->pageid) {
                            res.second = new key_type[index->valuesize];
                            key_type *inmemdata = inpage->inmemdata;
                            res.first = inmemdata[pageoff];
                            for (int ki = 0; ki < index->valuesize; ki++) {
                                res.second[ki] = inmemdata[pageoff + 1 + ki];
                            }
                            auto recordnum = pageoff/diskfilm->recordsize;
                            this->unset_bit(recordnum, index_pageid);
                        } else {
                            res = diskfilm->odirectreadfromdisk(
                                    pageid, pageoff,this, index_pageid);    // if readfromdisk indicating doesn't use o_direct;
                        }

                        index_res.findleaf->slotdata[index_res.slot] = index_res.findleaf->intrachain.put(
                                index_res.slot, res.second);
                        index_res.findleaf->locbitmap[index_res.slot] = true;
                        lru->put(index_res.findleaf->startkey, index_res.findleaf);// 更新 interchain
                        index->inkeynum++;
                        index->exkeynum--;
                    }
                }
                return res;
            }

        }

//        template< typename stat_type>
        bool delete_key(key_type key, disk_type* diskfilm){
//            struct timeval tqt1, tqt2, xt1, xt2, plt1, plt2, prdt1, prdt2;
//            double pltimeuse, prdtimeuse, pwdtimeuse, timeuse, xtimeuse;  // 读磁盘的时间，写磁盘的时间

//            gettimeofday(&xt1, NULL);


            auto [success_deleted,in_mem_flag, short_reduced_evict, m_size, del_piece] = index->delete_one(key);
//            query_stats->memnum += 1;
            if (success_deleted){
                if (in_mem_flag){
                    index->inkeynum--;
                    if (m_size == 0){ // 说明是在m_piece 中删除的数据, 并且intrachain size 为0
                        lru->remove(del_piece->startkey);
                    }
                }
                else{
                    index->exkeynum--;
//                query_stats->disknum += 1;
//                query_stats->diskpagenum += 1;
                    auto [pageid, offset,index_pageid] = query_flat_structure(short_reduced_evict,diskfilm->reduce_factor);
                    reduced_evicttable[index_pageid+1] -= 1;   // page中 的数据被删除掉一个，所以需要 the number of keys in page -1;
                    // 修改bitmap for records 表明, 该位置的record 已经invalid
                    auto recordnum = offset/diskfilm->recordsize;   // 计算是该页的第几个record
                    this->unset_bit(recordnum, index_pageid);
                }
                return true;
            }
            else{
                return false;
            }

            // 删除完成，以下的代码是为了测试，是否已经删除
            //            pair<key_type, key_type *> res;
//            auto index_res = index->search_one(key);
//            if (index_res.find == false) {   // 从 sort_piece 已经删除了数据，我们如何check: meiy delbitmap
//                if (index_res.flags) {  // 要删除的数据在内存中， 直接删除； 需要将b_piece 对应位置的key 修改为其后的 （next key）， 同时修改 bitmap_
//
//                    typename index_type::sortPiece* findpiece = (typename index_type::sortPiece*) (index_res.findleaf);
//                    auto finddata = (key_type *) (findpiece->ALEX_DATA_NODE_PAYLOAD_AT(index_res.slot));
//                    if (finddata != NULL)
//                        COUT_THIS("finddata != NULL ~ Lord, please come, chao needs You");
//                    if (findpiece->ALEX_DATA_NODE_KEY_AT(index_res.slot) == key)
//                        COUT_THIS("findpiece->ALEX_DATA_NODE_KEY_AT(index_res.slot) == key~~ "
//                                  "Lord, please come, chao needs You");
//                }
//                else {
//
//
//                    typename index_type::sortPiece* findpiece = (typename index_type::sortPiece*) (index_res.findleaf);
//                    if (findpiece->ALEX_DATA_NODE_KEY_AT(index_res.slot) == key)
//                        COUT_THIS("Lord, please come, chao needs You");
////                    pair<pageid_type, pageOff_type>  reduced_writeevict (0,0);
//
//                    auto short_reduced_evict = (unsigned long) findpiece->ALEX_DATA_NODE_PAYLOAD_AT(index_res.slot);  //原本这里可以直接是 writeevict 指向的是 pageid and offset
//                    // 在FILM+ flattened structures 中， 为进一步降低 memory usage，改成了flattened structure 的结构
////                    gettimeofday(&xt1, NULL);
//
//                    auto [pageid, offset,index_pageid] = query_flat_structure(short_reduced_evict,diskfilm->reduce_factor);
//
//                    // 将 querykey merge到memory，修改 bitmap，修改num in page
//
//                    reduced_evicttable[index_pageid+1] -= 1;
//
////                    gettimeofday(&prdt1, NULL);
//                    if (pageid == inpage->pageid) {
//                        res.second = new key_type[index->valuesize];
//                        key_type *inmemdata = inpage->inmemdata;
//                        res.first = inmemdata[offset];
//                        for (int ki = 0; ki < index->valuesize; ki++) {
//                            res.second[ki] = inmemdata[offset + 1 + ki];
//                        }
//                    } else {
//                        res = diskfilm->odirectreadfromdisk(pageid,offset,this,index_pageid);    // if readfromdisk indicating doesn't use o_direct;
//                    }
//
//                    findpiece->ALEX_DATA_NODE_PAYLOAD_AT(index_res.slot) = res.second;
//
//                    findpiece->ALEX_DATA_NODE_FLAG_AT(index_res.slot) = true;
//
//                    lru->put(findpiece->ALEX_DATA_NODE_KEY_AT(0), findpiece);// 更新 interchain
//
////                    gettimeofday(&plt2, NULL);
////                    pltimeuse =
////                            (plt2.tv_sec - plt1.tv_sec) + (double) (plt2.tv_usec - plt1.tv_usec) / 1000000.0;
////                    query_stats->lrutime += pltimeuse;
//
//                    index->inkeynum++;
//                    index->exkeynum--;
//                }
//            }
//            else {
//                if (index_res.flags){// 读取内存中已经删除的数据
//
//                    if (index_res.findleaf->delbitmap[index_res.slot] == false)
//                        COUT_THIS("Lord, please come, chao needs You")
//                }
//                else { //读取从磁盘删除的数据
//                    if (index_res.findleaf->delbitmap[index_res.slot] == false)
//                        COUT_THIS("Lord: I will remember the covenant I made with you in the days of your youth, "
//                                  "and I will establish an everlasting covenant with you!")
//                }
//            }
////            index_res = index->search_one(key);
//            return res;
        }

        template< typename stat_type>
        pair<key_type, key_type *> update_key(key_type key, key_type* new_payload, stat_type query_stats, disk_type diskfilm){
            struct timeval tqt1, tqt2, xt1, xt2, plt1, plt2, prdt1, prdt2;
            double pltimeuse, prdtimeuse, pwdtimeuse, timeuse, xtimeuse;  // 读磁盘的时间，写磁盘的时间
            pair<key_type, key_type *> res;
            gettimeofday(&xt1, NULL);
            auto index_res = index->search_one(key);
            gettimeofday(&xt2, NULL);
            xtimeuse = (xt2.tv_sec - xt1.tv_sec) + (double) (xt2.tv_usec - xt1.tv_usec) / 1000000.0;
            query_stats->xlookuptime += xtimeuse;
            if (index_res.find == false) {   // 从 sort_piece 中 读取数据
                if (index_res.flags) {
                    query_stats->memnum += 1;
                    typename index_type::sortPiece* findpiece = (typename index_type::sortPiece*) (index_res.findleaf);
                    findpiece->slotdata[index_res.slot] = new_payload;
//                    gettimeofday(&plt1, NULL);

                    // 更新 interchain
                    this->lru->put(findpiece->slotkey[0], findpiece);
//                    gettimeofday(&plt2, NULL);
//                    pltimeuse =
//                            (plt2.tv_sec - plt1.tv_sec) + (double) (plt2.tv_usec - plt1.tv_usec) / 1000000.0;
//                    query_stats->lrutime += pltimeuse;
                } else {
                    query_stats->disknum += 1;
                    query_stats->diskpagenum += 1;

                    typename index_type::sortPiece* findpiece = (typename index_type::sortPiece*) (index_res.findleaf);

//                    pair<pageid_type, pageOff_type>  reduced_writeevict (0,0);

                    auto short_reduced_evict = (unsigned long) findpiece->slotdata[index_res.slot];  //writeevict 指向的是 pageid and offset
                    auto int_reduced_evict = short_reduced_evict/2.0;
                    auto fraction = modf(int_reduced_evict,&int_reduced_evict);
                    auto int_high_low = this->reduced_evicttable[int_reduced_evict];
                    pageOff_type pageoff = 0;
                    // 根据reduced_evict 算出（读出） page id，算出读出offset
                    if (fraction){
                        pageoff = int_high_low&0xFFFF;
//                                unsigned short hi = hhh>>16 ;// 高16位
//                                unsigned short lo = hhh&0xFFFF ;// 低16位
                    }
                    else
                        pageoff =  int_high_low>>16 ;   // 高16位
                    unsigned long index_pageid = floor(int_reduced_evict/diskfilm.reduce_factor)*diskfilm.reduce_factor;
                    pageid_type pageid = this->reduced_evicttable[index_pageid];
                    // 将 querykey merge到memory，修改 bitmap，修改num in page

                    reduced_evicttable[index_pageid+1] -= 1;

//                    gettimeofday(&prdt1, NULL);
                    if (pageid == inpage->pageid) {
                        res.second = new key_type[index->valuesize];
                        key_type *inmemdata = inpage->inmemdata;
                        res.first = inmemdata[pageoff];
                        for (int ki = 0; ki < index->valuesize; ki++) {
                            res.second[ki] = inmemdata[pageoff + 1 + ki];
                        }
                    } else {
                        res = diskfilm.odirectreadfromdisk(
                                pageid, pageoff,this,index_pageid);    // if readfromdisk indicating doesn't use o_direct;
                    }
//                    gettimeofday(&prdt2, NULL);
//                    prdtimeuse = (prdt2.tv_sec - prdt1.tv_sec) +
//                                 (double) (prdt2.tv_usec - prdt1.tv_usec) / 1000000.0;
//                    query_stats->rdisktime += prdtimeuse;
/*
                    if (res.first != key) {
                        cout << "not find, what's wrong? my Lord, i need You~~~~" << key << endl;
                    }
*/
//                    gettimeofday(&plt1, NULL);

                    findpiece->slotdata[index_res.slot] = new_payload;
                    findpiece->locbitmap[index_res.slot] = true;

                    this->lru->put(findpiece->slotkey[0], findpiece);// 更新 interchain

//                    gettimeofday(&plt2, NULL);
//                    pltimeuse =
//                            (plt2.tv_sec - plt1.tv_sec) + (double) (plt2.tv_usec - plt1.tv_usec) / 1000000.0;
//                    query_stats->lrutime += pltimeuse;

                    index->inkeynum++;
                    index->exkeynum--;
                }
            } else {
                if (index_res.flags) {// 从内存中读数据
                    query_stats->memnum += 1;
                    auto finddata = (adalru::Node<lruOff_type, key_type *> *) index_res.findleaf->slotdata[index_res.slot];

                    finddata->value = new_payload;
//                    gettimeofday(&plt1, NULL);
                    index_res.findleaf->intrachain.moveTohead(finddata);// update intrachain

                    // 更新 interchain
                    lru->put(index_res.findleaf->startkey, index_res.findleaf);
//                    gettimeofday(&plt2, NULL);
//                    pltimeuse =
//                            (plt2.tv_sec - plt1.tv_sec) + (double) (plt2.tv_usec - plt1.tv_usec) / 1000000.0;
//                    query_stats->lrutime += pltimeuse;
                } else { //从磁盘读数据
                    query_stats->disknum += 1;
                    query_stats->diskpagenum += 1;

                    pair<pageid_type, pageOff_type>  reduced_writeevict (0,0);
//                            auto writeevict = (pair<pageid_type, pageOff_type> *) index_res.findleaf->slotdata[index_res.slot];  //writeevict 指向的是 pageid and offset
                    auto short_reduced_evict = (unsigned long) index_res.findleaf->slotdata[index_res.slot];
                    auto int_reduced_evict = short_reduced_evict/2.0;
                    auto fraction = modf(int_reduced_evict,&int_reduced_evict);
                    auto int_high_low = this->reduced_evicttable[int_reduced_evict];
                    // 根据reduced_evict 算出（读出） page id，算出读出offset
                    unsigned long index_pageid = floor(int_reduced_evict/diskfilm.reduce_factor)*diskfilm.reduce_factor;
                    pageid_type pageid = this->reduced_evicttable[index_pageid];
                    pageOff_type pageoff = 0;
                    if (fraction){
                        pageoff = int_high_low&0xFFFF;
//                                unsigned short hi = hhh>>16 ;// 高16位
//                                unsigned short lo = hhh&0xFFFF ;// 低16位
                    }
                    else
                        pageoff =  int_high_low>>16 ;   // 高16位
                    // 将 querykey merge到memory，修改 bitmap，修改num in page
//                            auto numinpage = memoryfilm.reduced_evicttable[index_pageid+1];
                    this->reduced_evicttable[index_pageid+1] -= 1;
//                    gettimeofday(&prdt1, NULL);
                    if (reduced_writeevict.first == inpage->pageid) {
                        res.second = new key_type[index->valuesize];
                        key_type *inmemdata = inpage->inmemdata;
                        res.first = inmemdata[reduced_writeevict.second];
                        for (int ki = 0; ki < index->valuesize; ki++) {
                            res.second[ki] = inmemdata[reduced_writeevict.second + 1 + ki];
                        }
                    } else {
                        res = diskfilm.odirectreadfromdisk(pageid, pageoff,this, index_pageid);    // if readfromdisk indicating doesn't use o_direct;
                    }
//                    gettimeofday(&prdt2, NULL);
//                    prdtimeuse = (prdt2.tv_sec - prdt1.tv_sec) +
//                                 (double) (prdt2.tv_usec - prdt1.tv_usec) / 1000000.0;
//                    query_stats->rdisktime += prdtimeuse;
/*
                    if (res.first != key) {
                        cout << "not find, what's wrong? my Lord, i need You~~~~" << key << endl;
                    }
                    */
//                    gettimeofday(&plt1, NULL);

                    index_res.findleaf->slotdata[index_res.slot] = index_res.findleaf->intrachain.put(
                            index_res.slot, new_payload);
                    index_res.findleaf->locbitmap[index_res.slot] = true;
                    lru->put(index_res.findleaf->startkey, index_res.findleaf);// 更新 interchain
//                    gettimeofday(&plt2, NULL);
//                    pltimeuse =
//                            (plt2.tv_sec - plt1.tv_sec) + (double) (plt2.tv_usec - plt1.tv_usec) / 1000000.0;
//                    query_stats->lrutime += pltimeuse;

                    index->inkeynum++;
                    index->exkeynum--;
                }
            }
            return res;
        }

        template< typename stat_type>
        pair<key_type, key_type *> scan_key(key_type key, int scan_num, stat_type query_stats, disk_type diskfilm){

        }

        template< typename stat_type>
        typename index_type::htaprange_res_find get_range(const key_type firstkey, const key_type lastkey, stat_type query_stats, disk_type diskfilm){
//            struct timeval tqt1, tqt2, xt1, xt2, plt1, plt2, prdt1, prdt2;
//            double pltimeuse, prdtimeuse, pwdtimeuse, timeuse, xtimeuse;  // 读磁盘的时间，写磁盘的时间
//            gettimeofday(&xt1, NULL);
            auto index_res = index->search_range(firstkey,lastkey);
//            gettimeofday(&xt2, NULL);
            // 根据 index_res 读取数据
            return index_res;


        }

        template< typename stat_type>
        typename index_type::htaprange_res_find get_scan(const key_type firstkey, int scan_num, stat_type query_stats, disk_type diskfilm){
//            struct timeval tqt1, tqt2, xt1, xt2, plt1, plt2, prdt1, prdt2;
//            double pltimeuse, prdtimeuse, pwdtimeuse, timeuse, xtimeuse;  // 读磁盘的时间，写磁盘的时间
//            gettimeofday(&xt1, NULL);
            auto index_res = index->search_scan(firstkey,scan_num);
//            gettimeofday(&xt2, NULL);
            // 根据 index_res 读取数据
            return index_res;


        }


        struct memoryusage
        {
            /// the total used memory
            double  totalusemem;
            /// the usage of each component
            infomap meminfo;
            /// Constructor
            inline memoryusage(double total,infomap eachmem)
                    : totalusemem(total), meminfo(eachmem)
            { }
            inline memoryusage()
                    : totalusemem(), meminfo()
            { }

        };

        struct inmempage  // a buffer page in memory
        {
            pageid_type pageid;
            //vector<key_type> inmemdata;
            key_type* inmemdata = NULL;
            int freespace;
            int usedspace;
            int recordnum;
            int record_num_reserved;
            inmempage(pageid_type idinpage,int sizepage,int numrecord){
                pageid = idinpage;
                freespace = sizepage;
                inmemdata = new key_type[sizepage];
                usedspace = 0;
                recordnum = 0;
                record_num_reserved = numrecord;
            }
            inmempage(){
            }
            ~inmempage(){
                delete [] inmemdata;
            }


        };

        // define the current page in memory, when the page is full, write to disk
        inmempage *inpage = NULL;
        vector<inmempage*> inmempages ;
        // 关于bitmap值的计算。 numrecord/8/4


        forceinline void createinmempage(int sizepage,int numrecord){
//            if (inpageid == 15437){
//                cout << "i need You, my Lord!! reused_pageid == 3027" << endl;
//            }
            inpage = new inmempage(inpageid,sizepage,numrecord);  //create a page in memory;
            inpageid += 1;
        }

        forceinline inmempage* createinmempage(unsigned int reused_pageid, int sizepage,int numrecord){
//            if (reused_pageid == 62002 ){
//                cout << "i need You, my Lord!!  reused_pageid == 3027" << endl;
//            }
            inmempage* reused_inpage = new inmempage(reused_pageid,sizepage,numrecord);  //create a page in memory;
            return reused_inpage;
        }

        // compute the usage of memory, that the total used memory by index,lru,data.addusage;


//        memoryusage runtimecomputeMemUsage(){
//
//            infomap indexdatausage = index->runtimeshow_verify();
//            double lruusage = indexdatausage["lruusage"] ;
//            double indexusage = indexdatausage["indexusage"] ;
//            double datausage = indexdatausage["datausage"];
//            double addusage = indexdatausage["addusage"] + double(reduced_evicttable.size()*sizeof(pageid_type))/1048576;
////            double newtracking =  (inpageid + 17 + index->exkeynum/2)*sizeof(pageid_type);
////            double newtracking1 = reduced_evicttable.size()*sizeof(pageid_type);
//            double totalmem = lruusage + datausage + indexusage + addusage;
//            double leaves = indexdatausage["leaves"];
//            double levels = indexdatausage["levels"];
//            double innernodes = indexdatausage["inners"];
//            double exkeynum = index->exkeynum;
//            double inkeynum = index->inkeynum;
//
//            infomap devied_mem;
//            devied_mem.insert(pair<std::string,double>("datausage",datausage));
//            devied_mem.insert(pair<std::string,double>("indexusage",indexusage));
//            devied_mem.insert(pair<std::string,double>("lruusage",lruusage));
//            devied_mem.insert(pair<std::string,double>("addusage",addusage));  //key,flag,pageID,offset
//            devied_mem.insert(pair<std::string,double>("leaves",leaves));
//            devied_mem.insert(pair<std::string,double>("levels",levels));
//            devied_mem.insert(pair<std::string,double>("innernodes",innernodes));
//            devied_mem.insert(pair<std::string,double>("totalusage",totalmem));  //key,flag,pageID,offset
//            devied_mem.insert(pair<std::string,double>("exkeynum",exkeynum));
//            devied_mem.insert(pair<std::string,double>("inkeynum",inkeynum));
//            memoryusage res_memusage(totalmem,devied_mem);
//            return res_memusage;
//
//        }
//

        memoryusage computeMemUsage(){

            infomap indexdatausage = index->show_verify();
            double hashlrusize = 24 +28;   // hash (key,v, next); doubly-link chain（k(unsigned int)，,v, prev,next）  这个更准确
            double globallruusage = hashlrusize * lru->size/double(1048576);
            double lruusage = indexdatausage["lruusage"]+globallruusage ;
            double indexusage = indexdatausage["indexusage"];
            double datausage = indexdatausage["datausage"];
            double pageusage =  double(reduced_evicttable.size()*sizeof(pageid_type))/1048576;
            double addusage = indexdatausage["addusage"] + pageusage;
            double preusage = indexdatausage["preusage"];
            double usepreusage = indexdatausage["usepreusage"];
            double totalmem = lruusage + datausage + indexusage + addusage +preusage-usepreusage;
            double leaves = indexdatausage["leaves"];
            double levels = indexdatausage["levels"];
            double innernodes = indexdatausage["inners"];
            double exkeynum = index->exkeynum;
            double inkeynum = index->inkeynum;

            infomap devied_mem;
            devied_mem.insert(pair<std::string,double>("datausage",datausage));
            devied_mem.insert(pair<std::string,double>("pageusage",pageusage));
            devied_mem.insert(pair<std::string,double>("preusage",preusage));
            devied_mem.insert(pair<std::string,double>("usepreusage",usepreusage));
            devied_mem.insert(pair<std::string,double>("indexusage",indexusage));
            devied_mem.insert(pair<std::string,double>("lruusage",lruusage));
            devied_mem.insert(pair<std::string,double>("addusage",addusage));  //key,flag,pageID,offset
            devied_mem.insert(pair<std::string,double>("leaves",leaves+index->newmodels));
            devied_mem.insert(pair<std::string,double>("buffers",leaves));
            devied_mem.insert(pair<std::string,double>("levels",levels));
            devied_mem.insert(pair<std::string,double>("innernodes",innernodes));
            devied_mem.insert(pair<std::string,double>("totalusage",totalmem));  //key,flag,pageID,offset
            devied_mem.insert(pair<std::string,double>("exkeynum",exkeynum));
            devied_mem.insert(pair<std::string,double>("inkeynum",inkeynum));
            memoryusage res_memusage(totalmem,devied_mem);
            return res_memusage;

        }

        pair<pageid_type ,pageOff_type>* writeevicttable(pair<pageid_type ,pageOff_type> pospage,pair<pageid_type ,pageOff_type>* evictpos){
            pair<pageid_type ,pageOff_type>* oldpospage = evictpos;
            evictpos->first = pospage.first;
            evictpos->second = pospage.second;
            return oldpospage;

//            evictpos->first = pospage.first;
//            evictpos->second = pospage.second;
        }

        pair<pageid_type ,pageOff_type>*  writeevicttable(pair<pageid_type ,pageOff_type> pospage){
            unsigned long int evictpos = evicttable.size();
//            evictkey.push_back(key);
            evicttable.push_back(pospage);
            auto eptr = &evicttable[evictpos];
            return eptr;
        }

        // Mark the entry for position in the bitmap
        forceinline void set_bit(int pos, unsigned long index_pageid){
            int bitmap_pos = pos >> 5;    // by chao, bitmap_pos represents the position of uint32
            int bit_pos = pos - (bitmap_pos << 5);     // by chao, bit_pos represents the position in uint32
            // 读出bitmap
            reduced_evicttable[index_pageid+2+bitmap_pos] |= (1U << bit_pos);
            // by chao, 将unsigned int 1 (32-bit 1), left shift bit_pos, 这表示，只将 bit_pos 位置为1， 其他位，采用 | 或操作，保持不变
//            vector<pageid_type> bitmap_int;
//            for (int i = 0; i<16;i ++){
//                bitmap_int.emplace_back(reduced_evicttable[index_pageid+2+i]);
//            }

//            bitmap_[bitmap_pos] |= (1ULL << bit_pos); // by chao, 将unsigned long long 1 (64-bit 1), left shift bit_pos, 这表示，只将 bit_pos 位置为1， 其他位，采用 | 或操作，保持不变
        }

        // Unmark the entry for position in the bitmap
        inline void unset_bit(int pos, unsigned long index_pageid) {
            int bitmap_pos = pos >> 5;
            int bit_pos = pos - (bitmap_pos << 5);
            reduced_evicttable[index_pageid+2+bitmap_pos]  &= ~(1ULL << bit_pos);

        }

        // Check whether the position corresponds to a key or has been merged into memory (as opposed to a hole)
        forceinline bool check_exists(int pos,  unsigned long index_pageid) const{
            int bitmap_pos = pos >> 5;  // divide by 32, by chao, since bitmap is a unsigned_int array
            int bit_pos = pos - (bitmap_pos << 5);
//            bool flag1 = static_cast<bool>(reduced_evicttable[index_pageid+2+bitmap_pos] & (1U << bit_pos));  // chao, need be removed s
            return static_cast<bool>(reduced_evicttable[index_pageid+2+bitmap_pos] & (1U << bit_pos));
        }

        unsigned long  reduced_writeevicttable(pageOff_type offset, int recordnum, int reduce_factor){
            unsigned long int evictpos = reduced_evicttable.size();
            // 判断一下是占用高字节还是低字节
            // 这里需要按照高/低字节位，以及short 类型对evictpos进行转换
            if ((recordnum & 1) == 0) { // 偶数, 占高字节
                // 将offset 放在高字节位置
                unsigned int iTest=0;
                unsigned short int *piTest=(unsigned short int *)&iTest;
//                *piTest=1024;	//低16位值
                piTest++;
                *piTest=offset;
                reduced_evicttable.push_back(iTest);
                // modify num in page for the current page
//                unsigned int num_in_page = reduced_evicttable[(evictpos/274)*274+1];
                unsigned long index_pageid = (evictpos/reduce_factor)*reduce_factor;
                reduced_evicttable[index_pageid+1] += 1;
                // modify bitmap for this page
                set_bit(recordnum,index_pageid);

                return (evictpos*2);
            }
            else{// 奇数，占低字节
                auto dual_off_iter = reduced_evicttable.end()-1 ;
                unsigned short int *piTest=(unsigned short int *)&(*dual_off_iter);
                *piTest=offset;	//低16位值
//                cout <<"the high 16 bits" << ((*dual_off_iter)>>16) << " the low 16 bits"<< ((*dual_off_iter)&0xFFFF)<< endl;
                // 或许可以利用inpage->pageid来计算

                // modify num in page for the current page
                evictpos -= 1;
                unsigned long index_pageid = ((evictpos)/reduce_factor)*reduce_factor;
                reduced_evicttable[index_pageid+1] += 1;
                // modify bitmap for this page
                set_bit(recordnum,index_pageid);

                return (evictpos*2+1);
            }


        }
        unsigned long  reduced_writeevicttable(pageOff_type offset, int recordnum, pageid_type pageid, int reduce_factor){
            unsigned long int evictpos = reduced_evicttable.size();
//            if (evictpos == 169822 || evictpos == 154383)
//                COUT_THIS("i need You, my Lord");
            // 判断一下是占用高字节还是低字节
            // 这里需要按照高/低字节位，以及short 类型对evictpos进行转换
            if ((recordnum & 1) == 0) { // 偶数, 占高字节
                // 将offset 放在高字节位置
                unsigned int iTest=0;
                unsigned short int *piTest=(unsigned short int *)&iTest;
//                *piTest=1024;	//低16位值
                piTest++;
                *piTest=offset;
//                auto xx = reduced_evicttable[154382];
                reduced_evicttable.push_back(iTest);

//                xx = reduced_evicttable[154382];
                // modify num in page for the current page
//                unsigned int num_in_page = reduced_evicttable[(evictpos/274)*274+1];
                unsigned long index_pageid = pageid*reduce_factor;
                reduced_evicttable[index_pageid+1] += 1;
                // modify bitmap for this page
                set_bit(recordnum,index_pageid);
//                xx = reduced_evicttable[154382];
                return (evictpos*2);
            }
            else{// 奇数，占低字节
                auto dual_off_iter = reduced_evicttable.end()-1 ;
                unsigned short int *piTest=(unsigned short int *)&(*dual_off_iter);
                *piTest=offset;	//低16位值
//                cout <<"the high 16 bits" << ((*dual_off_iter)>>16) << " the low 16 bits"<< ((*dual_off_iter)&0xFFFF)<< endl;
                // 或许可以利用inpage->pageid来计算

                // modify num in page for the current page
                evictpos -= 1;
                unsigned long index_pageid = pageid*reduce_factor;
//                auto xx = reduced_evicttable[154382];
                reduced_evicttable[index_pageid+1] += 1;
                // modify bitmap for this page
//                xx = reduced_evicttable[154382];
                set_bit(recordnum,index_pageid);
//                xx = reduced_evicttable[154382];
                return (evictpos*2+1);
            }


        }

        unsigned long  reduced_writeevicttable(pageid_type pageid, unsigned long index_pageid,pageOff_type offset, int recordnum,int reduce_factor){
            unsigned long int evictpos = pageid*reduce_factor+2+bitsnum+recordnum/2;
            // 判断一下是占用高字节还是低字节
            // 这里需要按照高/低字节位，以及short 类型对evictpos进行转换
            if ((recordnum & 1) == 0) { // 偶数, 占高字节
                // 将offset 放在高字节位置
                unsigned int iTest=0;
                unsigned short int *piTest=(unsigned short int *)&iTest;
                piTest++;
                *piTest=offset;
                // 修改对应位置的 高字节
                reduced_evicttable[evictpos] = iTest;
//                if ((pageid == 2798 && recordnum == 392) ) {
//                    cout <<"the high 16 bits" << ((iTest)>>16) << " the low 16 bits"<< ((iTest)&0xFFFF)<< endl;
//                    cout<< "i need You, my Lovely Lord" << endl; }
//                unsigned int num_in_page = reduced_evicttable[index_pageid+1];
//                unsigned long index_pageid_compute = (evictpos/reduce_factor)*reduce_factor;
                reduced_evicttable[index_pageid+1] += 1;
                // modify bitmap for this page
                set_bit(recordnum,index_pageid);

                return (evictpos*2);
            }
            else{// 奇数，占低字节
                // 修改对应位置的低字节
                auto dual_off_iter = reduced_evicttable.begin()+evictpos ;
                unsigned short int *piTest=(unsigned short int *)&(*dual_off_iter);
                *piTest=offset;	//低16位值

//                if ((pageid == 2798 && recordnum == 393) ) {
//                    cout <<"the high 16 bits " << ((*dual_off_iter)>>16) << " the low 16 bits "<< ((*dual_off_iter)&0xFFFF)<< endl;
//                    cout<< "i need You, my Lovely Lord" << endl;
//                }

                // modify num in page for the current page

//                unsigned long index_pageid = ((evictpos)/274)*274;
                reduced_evicttable[index_pageid+1] += 1;
                // modify bitmap for this page
                set_bit(recordnum,index_pageid);

                return (evictpos*2+1);
            }


        }

        // this function is used to check the num in a page, if the number is below a user-defined threshold
        // merge the page with the neighboring page？ (choose the second strategy) or merge the remaining keys into memory
        forceinline bool check_merge(unsigned long index_pageid){
//            auto num = reduced_evicttable[index_pageid+1];
            if (reduced_evicttable[index_pageid+1] < merge_threshold){
//                cout << "i need You, my Lord! merge the remaining keys" << endl;
                return true;
            }
//            cout << "thank You, my Lord" << endl;
            return false;
        }


        void page_merge(key_type *buf,unsigned int recordnum,unsigned long index_pageid){
            // 遍历bitmap, find the key still in page, and search the remaining key one by one
            // 记录page merge 的时间
//            struct timeval t1, t2;
//            double timeuse = 0.0;
//            gettimeofday(&t1, NULL);
//            int merge_num = 0;
            auto recordsize = this->index->valuesize+1;
            auto valuesize = this->index->valuesize;
            key_type value[valuesize];
            for (int i = 0; i < recordnum; i ++ ){   // (recordnum,int_reduced_evict,index_pageid
                if (this->check_exists(i, index_pageid)){ // if the bit is true, read the key from the read_page
                    // find the leaf piece that this key belonging to
                    key_type key = buf[i*recordsize];

                    for (int k = 0; k < valuesize; k ++)
                        value[k] = buf[i*recordsize+1+k];
//                    if (key == 9630175111229){
//                        COUT_THIS("thank You, my Lord, please be with my mom, key == 9630175111229");
//                    }
                    auto index_res = this->index->search_one(key);
                    // insert the key to the intrachain, 1. insert into the tail (choose this one), or 2. insert into the head
                    // 修改locbitmap         // 修改slotdata
                    if (index_res.findleaf->buffer_piece  != NULL){
//                        if (index_res.findleaf->slotkey[index_res.slot] != key)
//                            COUT_THIS("please come, my Lord!");
                        index_res.findleaf->slotdata[index_res.slot] = index_res.findleaf->intrachain.put(
                                index_res.slot, value);
                        index_res.findleaf->locbitmap[index_res.slot] = true;
                    }//cout << "thank You, my Lord! " << endl;
                    else{
                        typename index_type::sortPiece* findpiece = (typename index_type::sortPiece*) (index_res.findleaf);
//                        if ( findpiece->ALEX_DATA_NODE_KEY_AT(index_res.slot)!= key)
//                            COUT_THIS("please come, my Lord!");
                        findpiece->ALEX_DATA_NODE_PAYLOAD_AT(index_res.slot) = value;
                        findpiece->ALEX_DATA_NODE_FLAG_AT(index_res.slot) = true;
                    }
//                    merge_num ++; //  cout << "i need You, my Lord, merge this key" <<endl;
                }
            }
//            gettimeofday(&t2, NULL);
//            timeuse = (t2.tv_sec - t1.tv_sec) + (double) (t2.tv_usec - t1.tv_usec) / 1000000.0;
//            merge_timeuse += timeuse;
//            return merge_num;
        }

        void page_merge(key_type *buf,unsigned int recordnum,unsigned long index_pageid, unsigned long int offset){
            // 遍历bitmap, find the key still in page, and search the remaining key one by one
            // 记录page merge 的时间
            struct timeval t1, t2;
            double timeuse = 0.0;
            gettimeofday(&t1, NULL);
//            int merge_num = 0;
            auto recordsize = this->index->valuesize+1;
            auto valuesize = this->index->valuesize;
            key_type value[valuesize];
            for (int i = 0; i < recordnum; i ++ ){   // (recordnum,int_reduced_evict,index_pageid
                if (this->check_exists(i, index_pageid)){ // if the bit is true, read the key from the read_page
                    // find the leaf piece that this key belonging to
                    key_type key = buf[offset+i*recordsize];
//                    if ( key == 42266346608 )
//                        cout << "Jesus, where this key going?" << endl;
                    for (int k = 0; k < valuesize; k ++)
                        value[k] = buf[offset+i*recordsize+1+k];
                    auto index_res = this->index->search_one(key);
                    // insert the key to the intrachain, 1. insert into the tail (choose this one), or 2. insert into the head
                    // 修改locbitmap         // 修改slotdata
                    if (index_res.findleaf->buffer_piece  != NULL){
//                        if (index_res.findleaf->slotkey[index_res.slot] != key)
//                            COUT_THIS("please come, my Lord!");
                        index_res.findleaf->slotdata[index_res.slot] = index_res.findleaf->intrachain.put(
                                index_res.slot, value);
                        index_res.findleaf->locbitmap[index_res.slot] = true;
                    }//cout << "thank You, my Lord! " << endl;
                    else{
                        typename index_type::sortPiece* findpiece = (typename index_type::sortPiece*) (index_res.findleaf);

//                        if ( findpiece->ALEX_DATA_NODE_KEY_AT(index_res.slot)!= key)
//                            COUT_THIS("please come, my Lord!");
                        findpiece->ALEX_DATA_NODE_PAYLOAD_AT(index_res.slot) = value;
                        findpiece->ALEX_DATA_NODE_FLAG_AT(index_res.slot) = true;
                    }
//                    merge_num ++; //  cout << "i need You, my Lord, merge this key" <<endl;
                }
            }
            gettimeofday(&t2, NULL);
            timeuse = (t2.tv_sec - t1.tv_sec) + (double) (t2.tv_usec - t1.tv_usec) / 1000000.0;
            merge_timeuse += timeuse;
//            return merge_num;
        }

        void newpage_reduced_evicttable(pageid_type pageid){
//            if (pageid == 62002 ) {
//                auto index_pos = reduced_evicttable.size();
//                cout << "i need You, my Lord" << endl;
//            }
//            auto factor = reduced_evicttable.size()%138;
//
//            if (factor != 0)
//                COUT_THIS("Jesus, i need You");
            reduced_evicttable.push_back(pageid);
            // 写入num_in_page = 0
            reduced_evicttable.push_back(0);
            // 将全0 bitmap （8192个bool false，也就是256个 unsigned int 0），并写入bitmap
//            bitmap
            reduced_evicttable.insert(reduced_evicttable.end(),std::begin(bitmap),std::end(bitmap));
        }


        unsigned short int runtimeevictpagestodisk(disk_type *diskpage){
            unsigned short int pnum = inmempages.size();
            unsigned short int num = inmempages.size();  // 此次要写入磁盘的页的数目
            int fd = open(diskpage->pagefile, O_RDWR | O_DIRECT, 0755);
            unsigned long int fixed_buf_size = diskpage->pagesize * sizeof(key_type);  // 磁盘页固定的大小
            unsigned long seekdis =  fixed_buf_size * diskpage->nextpageid;
            lseek(fd, seekdis, SEEK_SET);
            int i = 0;
            while (pnum > diskpage->blocknum) {
                pnum -= diskpage->blocknum;    //diskpage->blocknum 表示 写出block 包含多少个page
                key_type *buf = new key_type[diskpage->pagesize* diskpage->blocknum];

                unsigned long buf_size = fixed_buf_size * diskpage->blocknum;
                int ret = posix_memalign((void **) &buf, 512, buf_size);

                int offset = 0;
                for(int k = 0;k<diskpage->blocknum;k++){
                    memcpy(buf + offset * diskpage->pagesize, &inmempages[i*diskpage->blocknum+k]->inmemdata[0], fixed_buf_size);
                    ++offset;
                }
                i += 1;
                ret = write(fd, buf, buf_size);

                delete [] buf;
            }

            key_type *buf = new key_type[diskpage->pagesize*pnum];
            unsigned long int buf_size = diskpage->pagesize * sizeof(key_type) * pnum;
            int ret = posix_memalign((void **) &buf, 512, buf_size);
            unsigned int offset = 0;
            for(int k = 0;k<pnum;k++){
                memcpy(buf + offset * diskpage->pagesize, &inmempages[i*diskpage->blocknum+k]->inmemdata[0], fixed_buf_size);
                ++offset;
            }
            ret = write(fd, buf, buf_size);
//            free(buf);
            delete[] buf;
            diskpage->nextpageid += num;
            close(fd);

            for (int mi = 0;mi < inmempages.size();mi++){
                delete inmempages[mi];
                inmempages[mi] = NULL;
            }
            inmempages.resize(0);
            return num;
        }

        unsigned short int runtimeevictreusedpagestodisk(disk_type *diskpage){
            unsigned short int num = inmempages.size();  // 此次要写入磁盘的页的数目
            int fd = open(diskpage->pagefile, O_RDWR | O_DIRECT, 0755);
            unsigned long int fixed_buf_size = diskpage->pagesize * sizeof(key_type);  // 磁盘页固定的大小
            auto initpageid = inmempages[0]->pageid;
            unsigned long seekdis = fixed_buf_size * initpageid;
            lseek(fd, seekdis, SEEK_SET);
            // lseek parameters, https://www.geeksforgeeks.org/lseek-in-c-to-read-the-alternate-nth-byte-and-write-it-in-another-file/
            // SEEK_SET it moves file pointer positon to the beginning of the file
            // SEEK_CUR specifies that the offset provided is relative to the current file position
            // SEEK_END it moves file pointer position to the end of file


            for (auto page : inmempages)
            {
                long off = long(page->pageid)- long(initpageid);
                long seekoff =  fixed_buf_size * off;
                lseek(fd, seekoff, SEEK_CUR);// the page id
                key_type *buf = new key_type[diskpage->pagesize];
                int ret = posix_memalign((void **) &buf, 512, fixed_buf_size);

                memcpy(buf , &page->inmemdata[0], fixed_buf_size);
                ret = write(fd, buf, fixed_buf_size);
                delete[] buf;
                initpageid = page->pageid+1;
            }

            close(fd);

            for (int mi = 0;mi < inmempages.size();mi++){
                delete inmempages[mi];
                inmempages[mi] = NULL;
            }
            vector<inmempage*>().swap(inmempages);
            return num;
        }

        pair<pageid_type ,pageOff_type>* evictkeytoinpage(key_type ekey,data_type edata, disk_type *diskpage){
//            if (ekey == 74431160)
//                cout << "Jesus, i need You" <<endl;
            if (!(inpage->recordnum --))  //如果还能容纳一个record，那么就向该页写出，如果不能，那么就将当前的inpage写出disk，创建新的inpage
            {
                inmempages.emplace_back(inpage);
//                diskpage->odirectenterpage(inpage->inmemdata);// write to disk
                createinmempage(diskpage->pagesize,diskpage->recordnum);
                inpage->recordnum -=1;

            }

            auto evictpos = writeevicttable(pair<pageid_type ,pageOff_type>(inpage->pageid, inpage->usedspace));
            inpage->inmemdata[inpage->usedspace++] = ekey;
            for (int vi = 0;vi < index->valuesize;vi++){
                inpage->inmemdata[inpage->usedspace++] = edata[vi];
            }

            return evictpos;

        }

        unsigned long reduced_evictkeytoinpage(key_type ekey,data_type edata, disk_type *diskpage){
//            if (ekey == 10289129908921){
//                cout << "Jesus, i need You, ekey == 18125561599" <<endl;
//            }

            if (!(inpage->freespace))  //如果还能容纳一个record，那么就向该页写出，如果不能，那么就将当前的inpage写出disk，创建新的inpage
            {
                inmempages.emplace_back(inpage);
                createinmempage(diskpage->pagesize,diskpage->recordnum);
                // 创建new page  时，需要对newpage 的information 进行初始化
                // 1, 将page id写如reduced_evicttable,  // 2, num in page 置为 0 // 3， bitmap 置0
                newpage_reduced_evicttable(inpage->pageid);
            }
            auto reduced_evictpos = reduced_writeevicttable(inpage->usedspace,inpage->recordnum,inpage->pageid,diskpage->reduce_factor);
            inpage->inmemdata[inpage->usedspace++] = ekey;
            inpage->recordnum +=1;
            for (int vi = 0;vi < index->valuesize;vi++){
                inpage->inmemdata[inpage->usedspace++] = edata[vi];
            }
            inpage->freespace -= diskpage->recordsize;
//            auto mmm = reduced_evicttable[8556414];
//            if (mmm != 0 && mmm != 62002)
//            {
//                COUT_THIS("thank You, my Lord!");
//            }
            return reduced_evictpos;
        }

        unsigned long reduced_evictkeytoinpage(key_type ekey,data_type edata, disk_type *diskpage, inmempage* reused_inpage, unsigned long index_pageid){
//            if (ekey == 10289129908921)
//                cout << "Jesus, i need You, ekey == 10289129908921 " <<endl;
            auto reduced_evictpos = reduced_writeevicttable(reused_inpage->pageid, index_pageid, reused_inpage->usedspace,reused_inpage->recordnum,diskpage->reduce_factor);
            reused_inpage->inmemdata[reused_inpage->usedspace++] = ekey;
            reused_inpage->recordnum +=1;
            for (int vi = 0;vi < index->valuesize;vi++){
                reused_inpage->inmemdata[reused_inpage->usedspace++] = edata[vi];
            }
            reused_inpage->freespace -= diskpage->recordsize;
            return reduced_evictpos;
        }


        void evictkeytoinpage(key_type ekey,data_type edata, disk_type *diskpage,pair<pageid_type ,pageOff_type>* evictpos){

            if (!(inpage->recordnum --))  //如果还能容纳一个record，那么就向该页写出，如果不能，那么就将当前的inpage写出disk，创建新的inpage
            {
                inmempages.emplace_back(inpage);
//                diskpage->odirectenterpage(inpage->inmemdata);// write to disk
                createinmempage(diskpage->pagesize,diskpage->recordnum);
                inpage->recordnum -=1;

            }
            writeevicttable(pair<pageid_type ,pageOff_type>(inpage->pageid, inpage->usedspace),evictpos);
            inpage->inmemdata[inpage->usedspace++] = ekey;
            for (int vi = 0;vi < index->valuesize;vi++){
                inpage->inmemdata[inpage->usedspace++] = edata[vi];
            }

        }



        unsigned long filmtransfer(unsigned int transleaves,disk_type *diskpage){   //perform the transfer procedure,
            //
            timeval initw1,initw2;
            unsigned long initwtime = 0;
            if (inpage == NULL){
                createinmempage(diskpage->pagesize,diskpage->recordnum);    //create a page in memory;
            }
            if (index->m_transleaf == NULL){
                index->m_transleaf = index->leaflevel.leafpieces[0];
            }
            /*
            if  ( this->index->Error > 256 && transleaves == 1 )
            {// 如果是这种情况，需要在一个leaf 中 批量写出，因为一次性写出leaf 所有的数据 会使得 内存的usage not full
                // 判断该页全部数据写出去 的usage
                auto midtransflag = simu_computeMemUsage(index->m_transleaf->slotkey.size());
                if (midtransflag.totalusemem > threshold ){
                    for (int k = 0; k < index->m_transleaf->slotkey.size();k++)
                    {
                        if (!(inpage->recordnum --))  //如果还能容纳一个record，那么就执行transfer，如果不能，那么就将当前的inpage写出disk，创建新的inpage
                        {
                            gettimeofday(&initw1, NULL);
                            diskpage->odirectenterpage(inpage->inmemdata);// write to disk , 这里是初始化阶段，所以没有采用 odirect 的方式 ，如果采用，只需要将在enterpage 前加上odirect， 去掉则不是直接操作磁盘
                            delete inpage;
                            inpage = NULL;
                            createinmempage(diskpage->pagesize,diskpage->recordnum);
                            gettimeofday(&initw2, NULL);
                            initwtime += (initw2.tv_sec - initw1.tv_sec) + (double) (initw2.tv_usec - initw1.tv_usec) / 1000000.0;
                            inpage->recordnum -=1;
                        }

                        auto eptr = writeevicttable(pair<unsigned long int,unsigned int>(inpage->pageid, inpage->usedspace));
                        inpage->inmemdata[inpage->usedspace++] = index->m_transleaf->slotkey[k];
                        auto transdata = (adalru::Node<unsigned int, key_type* >*) index->m_transleaf->slotdata[k];
                        for (int vi = 0; vi < index->valuesize; vi++){
                            inpage->inmemdata[inpage->usedspace++] = transdata->value[vi];
                        }

                        index->m_transleaf->intrachain.remove_node(transdata);//pop from intrachain
                        transdata = NULL;
                        index->m_transleaf->locbitmap[k] = false;
                        index->m_transleaf->slotdata[k]  = eptr ;  //
                    }
                    lru->remove(index->m_transleaf->startkey);// pop from interchain
                    index->inkeynum -= index->m_transleaf->slotkey.size();
                    index->exkeynum += index->m_transleaf->slotkey.size();
                    index->m_transleaf = index->leaflevel.leafpieces[index->leafsplit+1];//evict data from the first leaf
                    index->leafsplit += transleaves;
                }
                else{  // 一点一点的写出
                    auto midtrans = computeMemUsage();
                    int split = 0;
                    while (midtrans.totalusemem > threshold){
                        double ratio = (midtrans.totalusemem-threshold)/(midtrans.totalusemem-midtransflag.totalusemem ) ;
                        int trannum = ceil(index->m_transleaf->slotkey.size() * ratio+100);
                        for (int k = 0; k < trannum;k++)
                        {
                            if (!(inpage->recordnum --))  //如果还能容纳一个record，那么就执行transfer，如果不能，那么就将当前的inpage写出disk，创建新的inpage
                            {
                                gettimeofday(&initw1, NULL);
                                diskpage->odirectenterpage(inpage->inmemdata);// write to disk , 这里是初始化阶段，所以没有采用 odirect 的方式 ，如果采用，只需要将在enterpage 前加上odirect， 去掉则不是直接操作磁盘
                                delete inpage;
                                inpage = NULL;
                                createinmempage(diskpage->pagesize,diskpage->recordnum);
                                gettimeofday(&initw2, NULL);
                                initwtime += (initw2.tv_sec - initw1.tv_sec) + (double) (initw2.tv_usec - initw1.tv_usec) / 1000000.0;
                                inpage->recordnum -=1;
                            }

                            auto transdata = (adalru::Node<unsigned int, key_type* >*) index->m_transleaf->slotdata[split+k].second;
                            auto eptr = writeevicttable(pair<unsigned long int,unsigned int>(inpage->pageid, inpage->usedspace));
                            inpage->inmemdata[inpage->usedspace++] = index->m_transleaf->slotkey[split+k];

                            for (int vi = 0; vi < index->valuesize; vi++){
                                inpage->inmemdata[inpage->usedspace++] = transdata->value[vi];
                            }
                            index->m_transleaf->intrachain.remove_node(transdata);//pop from intrachain
                            transdata = NULL;
                            index->m_transleaf->slotdata[split+k].first = false;
                            index->m_transleaf->slotdata[split+k].second  = eptr ;  //
                        }
                        split += trannum;
                        index->inkeynum -= trannum;
                        index->exkeynum += trannum;
                        midtrans = computeMemUsage();
                    }
                }
            }
            else{
                */
            for (int i = 0; i< transleaves; i++){
                for (unsigned int k = 0; k < index->m_transleaf->slotkey.size();k++)
                {
                    if (!(inpage->recordnum --))  //如果还能容纳一个record，那么就执行transfer，如果不能，那么就将当前的inpage写出disk，创建新的inpage
                    {

                        gettimeofday(&initw1, NULL);
                        diskpage->odirectenterpage(inpage->inmemdata);// write to disk , 这里是初始化阶段，所以没有采用 odirect 的方式 ，如果采用，只需要将在enterpage 前加上odirect， 去掉则不是直接操作磁盘
                        gettimeofday(&initw2, NULL);
                        initwtime += (initw2.tv_sec - initw1.tv_sec) * 1000000 + (initw2.tv_usec - initw1.tv_usec);
                        delete inpage;
                        inpage = NULL;
                        createinmempage(diskpage->pagesize,diskpage->recordnum);
                        inpage->recordnum -=1;
                    }

                    auto eptr = writeevicttable(pair<pageid_type ,pageOff_type>(inpage->pageid, inpage->usedspace));
                    inpage->inmemdata[inpage->usedspace++] = index->m_transleaf->slotkey[k];
                    auto transdata = (adalru::Node<lruOff_type , key_type* >*) index->m_transleaf->slotdata[k];
                    for (int vi = 0; vi < index->valuesize; vi++){
                        inpage->inmemdata[inpage->usedspace++] = transdata->value[vi];
                    }

                    index->m_transleaf->intrachain.remove_node(transdata);//pop from intrachain
//                        delete transdata;
                    transdata = NULL;
                    index->m_transleaf->locbitmap[k] = false;
                    index->m_transleaf->slotdata[k]  = eptr ;  //
                }
                lru->remove(index->m_transleaf->startkey);// pop from interchain
                index->inkeynum -= index->m_transleaf->slotkey.size();
                index->exkeynum += index->m_transleaf->slotkey.size();
                index->m_transleaf = index->leaflevel.leafpieces[index->leafsplit+i+1];//evict data from the first leaf
            }
            index->leafsplit += transleaves;

            diskpage->initwtime += initwtime;
            return initwtime;
        }

        unsigned long filmtransfer(double totalusemem,disk_type *diskpage){   //perform the transfer procedure,
            //
            timeval initw1,initw2;
            unsigned long initwtime = 0;

            // 大致计算出此次要evict 的数量
            double ratio = (0.3) * (totalusemem - threshold) / threshold;
            unsigned int estimate_evict = ceil(index->inkeynum * ratio) + 300;
            unsigned int batch_evict = (1280 * diskpage->pagesize / (index->valuesize + 1));
            if (estimate_evict < batch_evict)
                batch_evict = estimate_evict;
            gettimeofday(&initw1, NULL);
            for (unsigned int i = 0; i < batch_evict; i++)   // the number of records to be evicted in batches
            {
                // evict data, get the tail node, remove the tail of intrachain
                auto evictleaf = lru->get_tail();
//                if (evictleaf->startkey == 139052311870037896){
//                    cout << "Jesus, wait for me!" << endl;
//                }
                auto evictslotV = evictleaf->intrachain.poptail();   // the tail of the accessed leaf
//                if (evictleaf->slotkey[evictslotV->key] == 142819719223789952){
//                    cout << "i need You, my lovely Lord, evicted_key == 142819719223789952" <<endl;
//                }
                auto writeevict = reduced_evictkeytoinpage(evictleaf->slotkey[evictslotV->key], evictslotV->value,
                                                           diskpage);
                evictleaf->locbitmap[evictslotV->key] = false;
                evictleaf->slotdata[evictslotV->key] = (void*)writeevict;
                delete evictslotV;
                evictslotV = NULL;

            }

            runtimeevictpagestodisk(diskpage);
            gettimeofday(&initw2, NULL);
            initwtime += (initw2.tv_sec - initw1.tv_sec)* 1000000 + (initw2.tv_usec - initw1.tv_usec);

            index->inkeynum -= batch_evict;
            index->exkeynum += batch_evict;
            diskpage->initwtime += initwtime;
            return initwtime;
        }



//        template<typename stat_type>
//        pair<bool,double > runtimejudgetrans(stat_type query_stats){
//            struct timeval ct1, ct2;
//            double ctimeuse;
//            gettimeofday(&ct1, NULL);
//            memoryusage res_memusage = runtimecomputeMemUsage();
//            gettimeofday(&ct2, NULL);
//            ctimeuse = (ct2.tv_sec - ct1.tv_sec) + (double) (ct2.tv_usec - ct1.tv_usec) / 1000000.0;
//            query_stats->computetimeuse += ctimeuse;
//            if (res_memusage.totalusemem > threshold) //perform transfer process
//                // 判断如果全部数据写出去都不能满足larger-than-memory data set
//                return pair<bool,double > (true,res_memusage.totalusemem);
//            else
//                return  pair<bool,double > (false,res_memusage.totalusemem);
//        }

        template<typename stat_type>
        pair<bool,memoryusage> judgetransfer(stat_type range_stats ){

            struct timeval ct1, ct2;
            double ctimeuse;
            gettimeofday(&ct1, NULL);
            memoryusage res_memusage = computeMemUsage();
            gettimeofday(&ct2, NULL);
            ctimeuse = (ct2.tv_sec - ct1.tv_sec) + (double) (ct2.tv_usec - ct1.tv_usec) / 1000000.0;

            map<string,double>::iterator iter;

//            for(iter = res_memusage.meminfo.begin(); iter != res_memusage.meminfo.end(); iter++)
//                cout<<iter->first<<" "<<iter->second<<" *** ";
//            cout<<endl;
            range_stats->computetimeuse += ctimeuse;
            if (res_memusage.totalusemem > threshold){//perform transfer process
                return {true,res_memusage};
            }
            return {false,res_memusage};
        }


        pair<bool,memoryusage> judgetransfer( ){
            memoryusage res_memusage;
            res_memusage = computeMemUsage();

            map<string,double>::iterator iter;

            for(iter = res_memusage.meminfo.begin(); iter != res_memusage.meminfo.end(); iter++)
                cout<<iter->first<<" "<<iter->second<<" *** ";
            cout<<endl;
            if (res_memusage.totalusemem > threshold){//perform transfer process
                return pair<bool,memoryusage>(true,res_memusage);
            }
            return  pair<bool,memoryusage>(false,res_memusage);
        }
    };

    typedef std::map<std::string,double> infomap;
    template <class key_type>
    class filmdisk{
    public:
        const char *pagefile;
        int pagesize;
        pageid_type nextpageid;
        int recordnum;  // 一个页最多容纳的 record  的数量
        int recordsize;
        int blocknum;
        int reduce_factor;
        unsigned long initwtime;
        filmdisk( const char* diskfile,int sizepage,int numrecord,int sizerecord, int merge_setting) {
            pagefile = diskfile;
            pagesize = sizepage;
            recordnum = numrecord;
            blocknum = 1024*1024/(sizepage*8);
            recordsize = sizerecord;
            reduce_factor = 2+(pagesize/recordsize)*17/32;
            nextpageid = 0;
            initwtime = 0.0;
        }


        vector<key_type> readfromdisk(pair<pageid_type ,pageOff_type> diskpos ,int sizerecord){   // use the buffer of memory
            FILE *fdisk;// 读取磁盘文件
            fdisk = fopen(pagefile,"rb+");
            fseek(fdisk,static_cast<unsigned long>(diskpos.first*pagesize*8),SEEK_SET);
            vector<key_type> pagedata;
            key_type tmp;

            for (int i = 0; i < pagesize; i ++){
                fread(&tmp, sizeof(long int), 1, fdisk); // 从文件中读数
                pagedata.push_back(tmp);
//                cout<< tmp << " ";
            }

            vector<key_type> res;

            for (int i = 0; i < sizerecord; i ++){
                res.push_back(pagedata[diskpos.second+i]);
            }

//            cout<<"Jesus, You are my refuge! "<<endl;
            fclose(fdisk);
            return res;

        }

        // use the system o_direct mode to read file
        vector<key_type> odirectreadfromdisk(pageid_type pageid ,pageOff_type offset ,int sizerecord){
            int fd;
            key_type *buf;
            vector<key_type> res;
            unsigned long int buf_size = pagesize*8;
            int ret = posix_memalign((void **)&buf, 512, buf_size);
            memset(buf, 'c', buf_size);

            fd = open(pagefile, O_RDWR | O_DIRECT , 0755);
            /*
            if (fd < 0){
                perror("open ./direct_io.data failed");
                exit(1);
            }
             */
            unsigned long seekdis = buf_size *pageid;
            ret = pread(fd, buf, buf_size,seekdis);
            if (ret <= 0){
                cout << "Jesus, i need You! ?????" << endl;
            }
            for (int i = 0; i < sizerecord; i ++){
                res.push_back(buf[offset+i]);
            }
//            cout<<"Jesus, You are my refuge! "<<endl;
//            free(buf);
            delete[] buf;
            close(fd);
            return res;
        }

        // use the system o_direct mode to read file
        pair<key_type,key_type*> odirectreadfromdisk(pair<pageid_type ,pageOff_type>* diskpos){
            int fd;
            key_type *buf;
            pair<key_type,key_type*> res;
            uint64_t buf_size = pagesize*8;
            int ret = posix_memalign((void **)&buf, 512, buf_size);
            memset(buf, 'c', buf_size);

            fd = open(pagefile, O_RDWR | O_DIRECT , 0755);
            uint64_t aoff = diskpos->first*buf_size;
            ret = pread(fd, buf, buf_size,aoff);
            if (ret <= 0){
                cout << "Jesus, i need You! ret <= 0" << endl;
            }

            res.first = buf[diskpos->second];
            res.second = new key_type[recordsize-1];
            for (int i = 0; i < recordsize-1; i ++){
                res.second[i] = buf[diskpos->second+1+i];
            }

            delete[] buf;
            close(fd);
            return res;
        }

        // use the system o_direct mode to read file
        template<class memory_type>
        pair<key_type,key_type*> odirectreadfromdisk(pageid_type pageid, pageOff_type offset, memory_type memoryfilm,
                                                     unsigned long index_pageid){
            int fd;
            key_type *buf;
            pair<key_type,key_type*> res;
            uint64_t buf_size = pagesize*8;
            int ret = posix_memalign((void **)&buf, 512, buf_size);
            memset(buf, 'c', buf_size);

            fd = open(pagefile, O_RDWR | O_DIRECT , 0755);
            uint64_t aoff = pageid *buf_size;
            ret = pread(fd, buf, buf_size,aoff);
            if (ret <= 0){
                cout << "Jesus, i need You! odirect read from disk = 0" << endl;
            }

            res.first = buf[offset];
            res.second = new key_type[recordsize-1];
            for (int i = 0; i < recordsize-1; i ++){
                res.second[i] = buf[offset+1+i];
            }

            // 如何确定，是page中的第几个key？
            auto recordnum = offset/this->recordsize;
//            bool former_flag = memoryfilm->check_exists(recordnum,index_pageid);
//            if (former_flag == false)
//                cout << "i need You, my Lord, the flag is wrong" << endl;

            /* when the query key is read from disk, need to unset the bitmap in the information of the page
            unset, that is let the key's bit set to 0, means that the key in page is out of date
             modify bitmap
            */
            memoryfilm->unset_bit(recordnum, index_pageid);
//            bool after_flag = memoryfilm->check_exists(recordnum,index_pageid);   // int pos, unsigned long evictpos, unsigned long index_pageid
//            if (after_flag == true)
//                cout << "i need You, my Lord! please come! unset bit wrong " << endl;
            // check merge
            if (memoryfilm->check_merge(index_pageid) == true){
                // record the index_pageid, that the position in the reduced_evictedtable, which can be reused.
                memoryfilm->reused_pageids.emplace_back(pageid );
                memoryfilm->num_page_merge += 1;
                memoryfilm->reused_index_pageid.emplace_back(index_pageid);
                memoryfilm->page_merge(buf,this->recordnum, index_pageid);
                memoryfilm->index->exkeynum -= memoryfilm->reduced_evicttable[index_pageid+1];
                memoryfilm->index->inkeynum += memoryfilm->reduced_evicttable[index_pageid+1];
                memoryfilm->reduced_evicttable[index_pageid+1] = 0;
//                cout << "Jesus, the last day of this year" << endl;
            }
            delete[] buf;
            close(fd);
            return res;
        }



        int enterpage(vector<key_type> datas)   //将一个内存页，写入磁盘, doesn't consider the odirect
        {

            FILE *fdisk;// 读取磁盘文件
            fdisk = fopen(pagefile,"rb+");
            fseek(fdisk,static_cast<unsigned long>(nextpageid * (pagesize)*8),SEEK_SET);
            nextpageid += 1;
//            cout<<datas.size()<<endl;
            for(int i = 0 ; i < datas.size() ; i++) {
                fwrite(&datas[i], sizeof(key_type), 1, fdisk);  //就是把id里面的值读到f里面，每次读8个字节，一共读size次
//                cout << datas[i] <<" ";
            }

            fflush(fdisk);
            fclose(fdisk);
//            cout<< "*********************** Jesus, finished one transfer*************************"<<endl;

        }


        void odirectenterpage(key_type* datas)   //将一个内存页，写入磁盘,consider the odirect

        {

            int fd;
            key_type *buf = new key_type[pagesize];

            unsigned long int buf_size = pagesize*sizeof(key_type);
            int ret = posix_memalign((void **)&buf, 512, buf_size);
//            memcpy(buf, &datas[0], datas.size()*sizeof(key_type));
//            int ret = posix_memalign((void **)&buf, 512, buf_size);
//            memset(buf, 'c', 4096);
            memcpy(buf, &datas[0],buf_size);

            fd = open(pagefile, O_RDWR | O_DIRECT , 0755);
            /*
            if (fd < 0){
                perror("open ./direct_io.data failed");
                exit(1);
            }
             */
            unsigned long seekdis = buf_size * nextpageid;
            ret = pwrite(fd, buf, buf_size,seekdis);
            if (ret <= 0){
                cout << "Jesus, i need You!" << endl;
            }
            nextpageid += 1;

//            cout<<"Jesus, You are my refuge! odirect write "<<endl;

//            free(buf);
            delete[] buf;
            close(fd);

//            cout<< "*********************** Jesus, finished one transfer*************************"<<endl;
        }

    };
}




#endif //EXPERIMENTCC12_FILMADASTORAGE_H
