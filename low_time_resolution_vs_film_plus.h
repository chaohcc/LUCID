//
// Created by CCMa on 2022/12/5.
//

#ifndef PLUSFILM_FILM_PLUS_H
#define PLUSFILM_FILM_PLUS_H

#include <stdlib.h>
#include <queue>
#include <sys/types.h>
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
//typedef double key_type;   //if (filename== "astro_ra")

#define forceinline inline __attribute__((__always_inline__))
using namespace std;


namespace filminsert{
#define PGM_SUB_EPS(x, epsilon) ((x) <= (epsilon) ? 0 : ((x) - (epsilon)))
#define PGM_SUB_EPS2(x, epsilon,size) ((x) <= (epsilon) || (x) - (epsilon)>=(size) ? 0 :((x) - (epsilon)) )
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
        double lrutime = 0.0;
        double disktime = 0.0;
        double rdisktime = 0.0;
        double wdisktime = 0.0;
        double sort_piecetime = 0.0;
        unsigned int pagecross = 0;  // 磁盘页 seek 的跨度
        double zipffactor = 0.0;
        string qtype = "point";
        string workload = "hotspot";  // "zipf", "random", "hotspot","randomzipf"
        double insert_frac = 0.0;
        double out_of_order_frac = 0.0;

        void print_stats(){
            struct timeval t1;
            gettimeofday(&t1, NULL);
            cout<<"ID " <<double(t1.tv_sec) << " qtype "<< qtype << " zipffactor "<< zipffactor << " insert_frac "<< insert_frac << " out_of_order_frac "<< out_of_order_frac << " memnum " << memnum << " disknum " << disknum << " crossnum "<< crossnum << " diskpagenum ";
            cout<< diskpagenum<<" crosspagenum " << crosspagenum <<" pagecross " << pagecross<<" writenum " << writenum<< " querytimeuse " <<querytimeuse <<" lrutimeuse " <<lrutime <<" disktimeuse " << disktime <<" readdisktimeuse " <<rdisktime
                <<" writedisktimeuse " <<wdisktime << " indexloopuptimeuse " << xlookuptime << " sort_piecetime "  << sort_piecetime << " workload " << workload <<"\n";
            ofstream savefile;
            string performance_file = "/home/wamdm/chaohong/clionDir/updatefilm/result/totally_rebuild_reduced_adalru_film_performance.txt";
            savefile.open(performance_file,ios::app);
//            savefile << "_ ";
            savefile<<"ID " <<double(t1.tv_sec) << " qtype "<< qtype  << " zipffactor "<< zipffactor << " insert_frac "<< insert_frac
                    << " out_of_order_frac "<< out_of_order_frac  << " memnum " << memnum << " disknum " << disknum << " crossnum "<< crossnum
                    << " diskpagenum ";
            savefile<< diskpagenum<<" crosspagenum " << crosspagenum <<" pagecross " << pagecross<<" writenum " << writenum
                    <<" querytimeuse " <<querytimeuse  <<" lrutimeuse " <<lrutime <<" disktimeuse " <<disktime << " readdisktimeuse " <<rdisktime
                    <<" writedisktimeuse " <<wdisktime << " indexloopuptimeuse " << xlookuptime<< " sort_piecetime "  << sort_piecetime << " workload " << workload<<"\n";
            savefile << flush;
            savefile.close();

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

    template <typename key_type,typename value_type>
    class FILMinsert{
    public:

        struct Leafpiece;  // 声明leafpiece 结构体
        struct sortPiece:Leafpiece{

            vector<key_type> slotkey;
//        vector< pair<bool, adalru::Node<key_type, std::vector<key_type> >*  > >slotdata;
            vector< bool> slotflag;   //adalru::Node<key_type, std::vector<key_type> >*sss
            vector< void* > slotdata;   //adalru::Node<key_type, std::vector<key_type> >*sss
            adalru::localLRU <lruOff_type,key_type* >  intrachain;    //unsigned short int

            sortPiece() = default;

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
                // 遍历 intrachain，所有 大于 end 的 值，都需要加1
                intrachain.modify(end);
//                if (end < 0 && this->slotkey.size() >0)// 所有 在 end 之后的，都需要加1
                this->slotkey.insert(slotkey.begin()+end+1,key);
                this->slotdata.insert( slotdata.begin()+end+1,intrachain.put(end+1,payload));
                this->slotflag.insert(slotflag.begin()+end+1,true);
                if (++this->intercept == this->slope ){
                    // rebuild this leaf piece, merge the slotkey and buffer piece
                    if (end < 0 && this->slotkey.size() >1){
                        cout << "i need You, my Lord, please come, rebuild this piece, and change the startkey" << endl;
                        return pair<bool,bool> (true,true);
                    }
                    else{
//                        cout << "i need You, my Lord, please come, rebuild this piece" << endl;
                        return pair<bool,bool> (true,false);
                    }
                }
                else{
                    if (end < 0 && this->slotkey.size() > 1){
//                        cout<< "i need You, my Lord" << endl;
                        return pair<bool,bool> (false,true);
                    }// 这表示start key 发生了变化
                    else{
                        return pair<bool,bool> (false,false);
                    }

                }

            }

            forceinline pair<bool,bool>  insert(key_type key, void* payload,bool flag) {
                this->slotkey.emplace_back(key);
                this->slotdata.emplace_back( payload);
                this->slotflag.emplace_back(flag);
                ++this->intercept;
            }

            // rebuild process 过程中遇到了 insert key, and slot data
            forceinline bool insert(key_type key,bool flag, void* payload) {

                int begin = 0, end = slotkey.size() - 1, mid;  // find last Equal or small, 找到最后一个大于或等于 key 的
                while (begin <= end) {
                    mid = (end + begin) / 2;
                    if (this->slotkey[mid] >= key) {
                        end = mid - 1;
                    } else
                        begin = mid + 1;
                }
                // 遍历 intrachain，所有 大于 end 的 值，都需要加1
                intrachain.modify(end);

                if (flag) {
                    this->slotdata.insert(slotdata.begin() + end + 1, intrachain.put(end + 1, (key_type *) payload));
                } else
                    this->slotdata.insert(slotdata.begin() + end + 1, payload);

                this->slotkey.insert(slotkey.begin() + end + 1, key);
                this->slotflag.insert(slotflag.begin() + end + 1, flag);
                if (++this->intercept == this->slope) {
                    cout << "i need You, my Lord, please come, rebuild this piece during rebuiding model" << endl;// rebuild this leaf piece, merge the slotkey and buffer piece
                }
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
            inline pair<bool,unsigned long> find_upper_in_leaf(const key_type &key,const int error){
                pair<bool,unsigned long> is_in(false,0);

                auto pos = (*this)(key);

                 is_in.first = true;
                 is_in.second = pos;

                 return is_in;
            }

        };
        struct Innerpiece;  // 声明leafpiece 结构体
        struct Innerlevel;  // 声明 innerlevel 结构体，每个innerllevel 包含一个 vector<Innerpiece>,一个 opt， 一个 pos；
        struct Leaflevel;
        unsigned int totalnum;
        unsigned int retrain = 0;  // 统计model retrain 的次数, 当该leaf 的buffer piece 达到阈值时
        int newmodels = 0;  // 由于retrain带来的多出来的new models的个数。只有当retrain 时models 个数增加，才会导致leaves 个数变化
        unsigned rebuild_inner = 0; // 记录inner level trebuild 的次数
        unsigned int inkeynum=0;
        unsigned int exkeynum=0;
        int valuesize = 0;
        int Error ;
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
                vector<key_type>().swap(node->slotkey);
                vector<bool>().swap(node->slotflag);
                vector<void*>().swap(node->slotdata);
                // release intrachain
                while (node->buffer_piece != NULL){
                    node = node->buffer_piece;
                    vector<bool>().swap(node->slotflag);
                    vector<void*>().swap(node->slotdata);
                }
                node->buffer_piece = NULL;
                delete node->buffer_piece;
                node->intrachain.deletelru();
//                for (int i = 0; i < node->intrachain.size; i ++){
//                    delete node->intrachain[i];
//                }
//                adalru::localLRU <lruOff_type,key_type* >  intrachain;    //unsigned short int
                delete node;
                node = NULL;
//                data_node_allocator().deallocate(static_cast<data_node_type*>(node), 1);
            }
        }

        void delete_inner(Innerpiece* node) {

            if (node = NULL)
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
            vector<FILMinsert<key_type,value_type>::Leafpiece*>().swap(leaflevel.leafpieces);
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
                m_tailleaf->slotflag.emplace_back( true);
            }
            else
            {
                auto a = leaflevel.opt->get_segment(m_tailleaf->endkey); // 将生成的new leaf piece插入到leaflevel 中
                m_tailleaf->update(a);
                interchain->put(m_tailleaf->startkey,m_tailleaf);
                auto buffersize = ceil(m_tailleaf->slotkey.size()*0.429);
                m_tailleaf->buffer_piece->slope = buffersize;
                m_tailleaf->buffer_piece->slotkey.reserve(buffersize);
                m_tailleaf->buffer_piece->slotdata.reserve(buffersize);
                m_tailleaf->buffer_piece->slotflag.reserve(buffersize);
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
                m_tailleaf->buffer_piece = new sortPiece();
                buffersize = ceil(8192*2*0.429);
                m_tailleaf->buffer_piece->slope = buffersize;
                m_tailleaf->buffer_piece->slotkey.reserve(buffersize);
                m_tailleaf->buffer_piece->slotdata.reserve(buffersize);
                m_tailleaf->buffer_piece->slotflag.reserve(buffersize);
                m_tailleaf->slotkey.reserve(8192*2);
                leaflevel.opt->append_point(p.first, 0);
                m_tailleaf->slotdata.emplace_back( m_tailleaf->intrachain.put(0,payload));
                m_tailleaf->slotkey.push_back(p.first);
                m_tailleaf->slotflag.emplace_back( true);
                m_tailleaf->startkey = p.first;
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
            if (innerlevels.size()>1)
                root = innerlevels.back()->innerpieces[0];
            this->verify_film(&vstats);
            cout<<"index append test finished,Jesus, my Lord, i need You, praise to You for ever and ever!" << endl;
        }

        /*
        inline void update_random(std::vector<key_type> keys,std::vector<key_type> payload){
            // 首先search , 如果不存在，则插入，如果存在，则 pass
            for (unsigned int i = 0; i < keys.size(); i ++){
                this->update_search_one(keys[i],payload);
            }
            this->inkeynum += keys.size();
        }
        */

        template<class lru_type>
        inline void insert_random(std::vector<key_type> keys,key_type* payload,lru_type* interchain){
            // 首先search , 如果不存在，则插入，如果存在，则 pass
            unsigned int i = 0;
            if (keys[i] == 62129282194 ){
                cout<< "Jesus, please come!, 161142159158" <<endl;
            }
            for (; i < keys.size(); i ++){
                this->insert_one(keys[i],payload,interchain);
            }
            this->inkeynum += keys.size();

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
            unsigned int slotleaf0;
            unsigned int slotleaf1;

            unsigned int slotbuffer0;
            unsigned int slotbuffer1;

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

        struct AlexCompare  {

            template <class T1, class T2>
            bool operator()(const T1 x,const T2 y) const {
                static_assert(
                        std::is_arithmetic<T1>::value && std::is_arithmetic<T2>::value,
                        "Comparison types must be numeric.");
                return x < y;
            }


//            template <class T1, class T2>
//            AlexCompare(const T1 x, const T2 y)(){
//
//                    static_assert(
//                            std::is_arithmetic<T1>::value && std::is_arithmetic<T2>::value,
//                            "Comparison types must be numeric.");
//                    return (x < y);
//
//
//            };

        };

        // True if a < b
        template <class K>
        forceinline bool key_less(const key_type a, const K b) const {
            return (a<b);
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

        // True if a >= b
        template <class K>
        forceinline bool key_greaterequal( key_type a, K b) const {
            return !(a<b);
        }

        // True if a == b
        template <class K>
        forceinline bool key_equal(const key_type a, const K b) const {
            return !(a<b) && !(b<a);
        }


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
                    key_type mkey = (*(lo+m-bound))->startkey  ;// chao
//        key_greater(mkey, key);
                    auto mkeydiff = mkey - key;  //chao
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



        forceinline result_find search_one(const key_type key,access_stats *query_stats) {
            // 从 root 开始，即，从root开始的结尾开始
            struct timeval ct1, ct2;
            double ctimeuse;
            if (vstats.innernum > 1){
                auto itlevel = innerlevels.end() -1; // root level 的 iterator
                // internal level use exponential search, the parameters used in exponential_search are: predict_pos, search key
                /*
                 * 1. get the predicted pos
                 * 2. get the [l,r] used in binary search or, get the [l,r]
                 * 3. use the binary search in [l,r] to get the final position
                 */

                auto predicted_pos = (*root)(key);
                auto childpieces = (*(--itlevel))->innerpieces;
                auto search_pos = exponential_search_upper_bound(childpieces,predicted_pos, key) -1;
                // 判断 找到的inner node 是否满足条件
                auto curit = (*(itlevel))->innerpieces.begin() + search_pos;
                FILMinsert<key_type,value_type>::Innerpiece* b = *curit;
//                if (b->startkey > key )
//                    cout << "i need You, my lovely Lord, b->startkey > key " << endl;

                for (; --itlevel >= innerlevels.begin();){
                    predicted_pos = (*b)(key);
                    childpieces = (*(itlevel))->innerpieces;
                    search_pos = exponential_search_upper_bound(childpieces,predicted_pos, key) -1;
                    // 判断 找到的inner node 是否满足条件
                    curit = (*(itlevel))->innerpieces.begin() + search_pos;
                }

                // find the leaf level
//                auto curleaf = leaflevel.leafpieces.begin();

                b = *curit;
                predicted_pos = (*b)(key);

                search_pos = exponential_search_upper_bound(leaflevel.leafpieces,predicted_pos, key) -1;
                // 预测位置和最终位置的距离
//                auto distance = abs(long (search_pos)-long(predicted_pos));
//                if (distance >= 512){
//                    cout << "i need You, my Lord, please come~~" << distance << endl;
//                }
                auto curleaf = leaflevel.leafpieces.begin() + search_pos;

                Leafpiece* a = *curleaf;

                // 确定 该leaf 的哪一个链 leaf piece
                while (key > a->endkey && a->buffer_piece) {
                    a = a->buffer_piece;
                }
                if (a->buffer_piece == NULL) {
                    sortPiece *buffer = (sortPiece *) a;
                    auto sortpos = (*buffer)(key);
                    if (buffer->slotkey[sortpos] == key) {
                        result_find index_res = result_find(false, buffer->slotflag[sortpos], buffer,
                                                            sortpos);   // sort_list data
                        return index_res;
                    } else {
                        cout << " i need You, my Lord, please help me! not find the key " << key << endl;
                    }
                }
                else {
                    predicted_pos = (*a)(key);

                    int slotlo = PGM_SUB_EPS2(predicted_pos, ErrorRecursive + 1, (a->slotkey.size() - 1));
                    int slothi = PGM_ADD_EPS(predicted_pos, ErrorRecursive, (a->slotkey.size() - 1));
                    int resslot = slotlo;
                    for (; resslot < slothi; resslot++) {
                        if (key <= a->slotkey[resslot]) {
                            break;
                        }
                    }

                    gettimeofday(&ct1, NULL);
                    if (key != a->slotkey[resslot]) {
                        // find in sort_list
                        while (a->buffer_piece->buffer_piece != NULL) {
                            a = a->buffer_piece;
                        }
                        sortPiece *buffer = (sortPiece *) a->buffer_piece;
                        auto sortpos = (*buffer)(key);
                        if (buffer->slotkey[sortpos] == key) {
                            result_find index_res = result_find(false, buffer->slotflag[sortpos], buffer,
                                                                sortpos);   // sort_list data
                            return index_res;
                        } else {
                            cout << " i need You, my Lord, please help me! not find the key " << key << endl;
                        }
                    }
                    gettimeofday(&ct2, NULL);
                    ctimeuse = (ct2.tv_sec - ct1.tv_sec) + (double) (ct2.tv_usec - ct1.tv_usec) / 1000000.0;
                    query_stats->sort_piecetime += ctimeuse;
                    result_find index_res = result_find(true, a->slotflag[resslot], a, resslot);   // regular data
                    return index_res;
                }
            }
            else{

                // 只有一层leaf
                static constexpr size_t linear_search_threshold = 8 * 64 / sizeof(Innerpiece);   //sizeof(Innerpiece)

                // leaf level
                auto curleaf = leaflevel.leafpieces.begin();
                auto predicted_pos = (*root)(key);

                auto leaflo = PGM_SUB_EPS2(predicted_pos, Error + 1,leaflevel.leafpieces.size()) + leaflevel.leafpieces.begin();
                auto leafhi = PGM_ADD_EPS(predicted_pos, Error,leaflevel.leafpieces.size()-1) + leaflevel.leafpieces.begin();
                auto h = PGM_ADD_EPS(predicted_pos, Error,leaflevel.leafpieces.size()-1);
//
//                Leafpiece*  y= *leafhi;
//                Leafpiece*  z= *leaflo;
                if (key >= (*leafhi)->startkey){curleaf = leafhi;}
                else{
                    if constexpr (4 <= linear_search_threshold ){   //ErrorRecursive
                        for (; (*leaflo)->endkey < key; leaflo++){
//                            Leafpiece*  x= *leaflo;
                            continue;
                        }
                        curleaf = leaflo;
                    }
                    else{
                        curleaf = std::lower_bound(leaflo, leafhi, key);   // upper_bound 返回在前闭后开区间查找key 的上届，返回大于 key 的第一个元素位置； lower_bound 返回第一个大于或等于 key 的元素的位置
                    }
                }
                Leafpiece* a = *curleaf;
                // 确定 该leaf 的哪一个链 leaf piece
                while (key > a->endkey && a->buffer_piece) {
                    a = a->buffer_piece;
                }
                if (a->buffer_piece == NULL) {
                    sortPiece *buffer = (sortPiece *) a;
                    auto sortpos = (*buffer)(key);
                    if (buffer->slotkey[sortpos] == key) {
                        result_find index_res = result_find(false, buffer->slotflag[sortpos], buffer,
                                                            sortpos);   // sort_list data
                        return index_res;
                    } else {
                        cout << " i need You, my Lord, please help me! not find the key " << key << endl;
                    }
                }
                else {
                    predicted_pos = (*a)(key);

                    int slotlo = PGM_SUB_EPS2(predicted_pos, ErrorRecursive + 1, (a->slotkey.size() - 1));
                    int slothi = PGM_ADD_EPS(predicted_pos, ErrorRecursive, (a->slotkey.size() - 1));
                    int resslot = slotlo;
                    for (; resslot < slothi; resslot++) {
                        if (key <= a->slotkey[resslot]) {
                            break;
                        }
                    }

                    gettimeofday(&ct1, NULL);
                    if (key != a->slotkey[resslot]) {
                        // find in sort_list
                        while (a->buffer_piece->buffer_piece != NULL) {
                            a = a->buffer_piece;
                        }
                        sortPiece *buffer = (sortPiece *) a->buffer_piece;
                        auto sortpos = (*buffer)(key);
                        if (buffer->slotkey[sortpos] == key) {
                            result_find index_res = result_find(false, buffer->slotflag[sortpos], buffer,
                                                                sortpos);   // sort_list data
                            return index_res;
                        } else {
                            cout << " i need You, my Lord, please help me! not find the key " << key << endl;
                        }
                    }
                    gettimeofday(&ct2, NULL);
                    ctimeuse = (ct2.tv_sec - ct1.tv_sec) + (double) (ct2.tv_usec - ct1.tv_usec) / 1000000.0;
                    query_stats->sort_piecetime += ctimeuse;
                    result_find index_res = result_find(true, a->slotflag[resslot], a, resslot);   // regular data
                    return index_res;
                }
            }
        }

        forceinline range_key_result_find search_one_lower(key_type key) {
            // 找到 第一个大于或者等于 search_key 的piece (包括leaf，buffer)， 及对应的位置
            struct timeval ct1, ct2;
            double ctimeuse;
            if (vstats.innernum > 1) {
                auto itlevel = innerlevels.end() - 1; // root level 的 iterator
                auto predicted_pos = (*root)(key);
                auto childpieces = (*(--itlevel))->innerpieces;
                auto search_pos = exponential_search_upper_bound(childpieces, predicted_pos, key) - 1;
                // 判断 找到的inner node 是否满足条件
                auto curit = (*(itlevel))->innerpieces.begin() + search_pos;
                FILMinsert<key_type, value_type>::Innerpiece *b = *curit;
//                if (b->startkey > key )
//                    cout << "i need You, my lovely Lord, b->startkey > key " << endl;

                for (; --itlevel >= innerlevels.begin();) {
                    predicted_pos = (*b)(key);
                    childpieces = (*(itlevel))->innerpieces;
                    search_pos = exponential_search_upper_bound(childpieces, predicted_pos, key) - 1;
                    curit = (*(itlevel))->innerpieces.begin() + search_pos;
                }

                b = *curit;
                predicted_pos = (*b)(key);
                search_pos = exponential_search_upper_bound(leaflevel.leafpieces, predicted_pos, key) - 1;
                auto curleaf = leaflevel.leafpieces.begin() + search_pos;

                Leafpiece *a = *curleaf;

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
                    index_res.slots.push_back(leaf_a->slotkey.size() - 1);
                    index_res.findleafs.emplace_back(a);
                    sortPiece *buffer = (sortPiece *) a;
                    auto sortpos = (*buffer)(key);
                    index_res.slots.push_back(sortpos);
                    return index_res;
                }
                else {
                    predicted_pos = (*a)(key);

                    int slotlo = PGM_SUB_EPS2(predicted_pos, ErrorRecursive + 1, (a->slotkey.size() - 1));
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
                    sortPiece *buffer = (sortPiece *) a->buffer_piece;
                    auto sortpos = (*buffer)(key);
                    if (sortpos >0 && sortpos == buffer->intercept)
                        sortpos -= 1;

                    index_res.findleafs.emplace_back(a->buffer_piece);
                    index_res.slots.push_back(sortpos);
                    return index_res;
                    }
                }
                /*
                else{

                    // 只有一层leaf
                    static constexpr size_t linear_search_threshold = 8 * 64 / sizeof(Innerpiece);   //sizeof(Innerpiece)

                    // leaf level
                    auto curleaf = leaflevel.leafpieces.begin();
                    auto predicted_pos = (*root)(key);

                    auto leaflo = PGM_SUB_EPS2(predicted_pos, Error + 1,leaflevel.leafpieces.size()) + leaflevel.leafpieces.begin();
                    auto leafhi = PGM_ADD_EPS(predicted_pos, Error,leaflevel.leafpieces.size()-1) + leaflevel.leafpieces.begin();
                    auto h = PGM_ADD_EPS(predicted_pos, Error,leaflevel.leafpieces.size()-1);
    //
    //                Leafpiece*  y= *leafhi;
    //                Leafpiece*  z= *leaflo;
                    if (key >= (*leafhi)->startkey){curleaf = leafhi;}
                    else{
                        if constexpr (4 <= linear_search_threshold ){   //ErrorRecursive
                            for (; (*leaflo)->endkey < key; leaflo++){
    //                            Leafpiece*  x= *leaflo;
                                continue;
                            }
                            curleaf = leaflo;
                        }
                        else{
                            curleaf = std::lower_bound(leaflo, leafhi, key);   // upper_bound 返回在前闭后开区间查找key 的上届，返回大于 key 的第一个元素位置； lower_bound 返回第一个大于或等于 key 的元素的位置
                        }
                    }
                    Leafpiece* a = *curleaf;
                    // 确定 该leaf 的哪一个链 leaf piece
                    while (key > a->endkey && a->buffer_piece) {
                        a = a->buffer_piece;
                    }
                    if (a->buffer_piece == NULL) {
                        sortPiece *buffer = (sortPiece *) a;
                        auto sortpos = (*buffer)(key);
                        if (buffer->slotkey[sortpos] == key) {
                            result_find index_res = result_find(false, buffer->slotflag[sortpos], buffer,
                                                                sortpos);   // sort_list data
                            return index_res;
                        } else {
                            cout << " i need You, my Lord, please help me! not find the key " << key << endl;
                        }
                    }
                    else {
                        predicted_pos = (*a)(key);

                        int slotlo = PGM_SUB_EPS2(predicted_pos, ErrorRecursive + 1, (a->slotkey.size() - 1));
                        int slothi = PGM_ADD_EPS(predicted_pos, ErrorRecursive, (a->slotkey.size() - 1));
                        int resslot = slotlo;
                        for (; resslot < slothi; resslot++) {
                            if (key <= a->slotkey[resslot]) {
                                break;
                            }
                        }

                        gettimeofday(&ct1, NULL);
                        if (key != a->slotkey[resslot]) {
                            // find in sort_list
                            while (a->buffer_piece->buffer_piece != NULL) {
                                a = a->buffer_piece;
                            }
                            sortPiece *buffer = (sortPiece *) a->buffer_piece;
                            auto sortpos = (*buffer)(key);
                            if (buffer->slotkey[sortpos] == key) {
                                result_find index_res = result_find(false, buffer->slotflag[sortpos], buffer,
                                                                    sortpos);   // sort_list data
                                return index_res;
                            } else {
                                cout << " i need You, my Lord, please help me! not find the key " << key << endl;
                            }
                        }
                        gettimeofday(&ct2, NULL);
                        ctimeuse = (ct2.tv_sec - ct1.tv_sec) + (double) (ct2.tv_usec - ct1.tv_usec) / 1000000.0;
                        query_stats->sort_piecetime += ctimeuse;
                        result_find index_res = result_find(true, a->slotflag[resslot], a, resslot);   // regular data
                        return index_res;
                    }
                }
                */
            }

        forceinline result_find search_one_leaf(const key_type key) {
            // 从 root 开始，即，从root开始的结尾开始
            struct timeval ct1, ct2;
            double ctimeuse;
            if (vstats.innernum > 1){
                auto itlevel = innerlevels.end() -1; // root level 的 iterator
                // internal level use exponential search, the parameters used in exponential_search are: predict_pos, search key
                /*
                 * 1. get the predicted pos
                 * 2. get the [l,r] used in binary search or, get the [l,r]
                 * 3. use the binary search in [l,r] to get the final position
                 */

                auto predicted_pos = (*root)(key);
                auto childpieces = (*(--itlevel))->innerpieces;
                auto search_pos = exponential_search_upper_bound(childpieces,predicted_pos, key) -1;
                // 判断 找到的inner node 是否满足条件
                auto curit = (*(itlevel))->innerpieces.begin() + search_pos;
                FILMinsert<key_type,value_type>::Innerpiece* b = *curit;
//                if (b->startkey > key )
//                    cout << "i need You, my lovely Lord, b->startkey > key " << endl;

                for (; --itlevel >= innerlevels.begin();){
                    predicted_pos = (*b)(key);
                    childpieces = (*(itlevel))->innerpieces;
                    search_pos = exponential_search_upper_bound(childpieces,predicted_pos, key) -1;
                    // 判断 找到的inner node 是否满足条件
                    curit = (*(itlevel))->innerpieces.begin() + search_pos;
                }

                // find the leaf level
//                auto curleaf = leaflevel.leafpieces.begin();

                b = *curit;
                predicted_pos = (*b)(key);

                search_pos = exponential_search_upper_bound(leaflevel.leafpieces,predicted_pos, key) -1;
                auto curleaf = leaflevel.leafpieces.begin() + search_pos;

                Leafpiece* a = *curleaf;

                result_find index_res= result_find(false,a->slotflag[0],a,0);   // not find, then, insert into bufferpiece,
                index_res.leafslot = search_pos;

                return index_res;
            }
            else{

                // 只有一层leaf
                static constexpr size_t linear_search_threshold = 8 * 64 / sizeof(Innerpiece);   //sizeof(Innerpiece)

                // leaf level
                auto curleaf = leaflevel.leafpieces.begin();
                auto pos = (*root)(key);

                auto leaflo = PGM_SUB_EPS2(pos, Error + 1,leaflevel.leafpieces.size()) + leaflevel.leafpieces.begin();
                auto leafhi = PGM_ADD_EPS(pos, Error,leaflevel.leafpieces.size()-1) + leaflevel.leafpieces.begin();
                auto h = PGM_ADD_EPS(pos, Error,leaflevel.leafpieces.size()-1);
//
//                Leafpiece*  y= *leafhi;
//                Leafpiece*  z= *leaflo;
                if (key >= (*leafhi)->startkey){curleaf = leafhi;}
                else{
                    if constexpr (4 <= linear_search_threshold ){   //ErrorRecursive
                        for (; (*leaflo)->endkey < key; leaflo++){
//                            Leafpiece*  x= *leaflo;
                            continue;
                        }

                        curleaf = leaflo;
                    }
                    else{ curleaf = std::lower_bound(leaflo, leafhi, key);   // upper_bound 返回在前闭后开区间查找key 的上届，返回大于 key 的第一个元素位置； lower_bound 返回第一个大于或等于 key 的元素的位置
                    }
                }
                Leafpiece* a = *curleaf;
                pos = (*a)(key);

                int slotlo = PGM_SUB_EPS(pos, ErrorRecursive + 1) ;  //,(--itlevel)->size()
                int slothi = PGM_ADD_EPS(pos, ErrorRecursive,(a->slotkey.size()-1));
                int resslot = slotlo;
                for (; resslot<slothi;resslot++){
                    if (key <= a->slotkey[resslot])
                    {
                        break;
                    }
                }
                result_find index_res;
                if (key != a->slotkey[resslot]){
                    index_res = result_find(false,a->slotflag[resslot],a,resslot);   // regular data
                }
                else
                    index_res = result_find(true,a->slotflag[resslot],a,resslot);   // regular data
                return index_res;
            }
        }

        forceinline result_find search_one(const key_type key) {
            // 从 root 开始，即，从root开始的结尾开始
            struct timeval ct1, ct2;
            double ctimeuse;
            if (vstats.innernum > 1){
                auto itlevel = innerlevels.end() -1; // root level 的 iterator
                // internal level use exponential search, the parameters used in exponential_search are: predict_pos, search key
                /*
                 * 1. get the predicted pos
                 * 2. get the [l,r] used in binary search or, get the [l,r]
                 * 3. use the binary search in [l,r] to get the final position
                 */

                auto predicted_pos = (*root)(key);
                auto childpieces = (*(--itlevel))->innerpieces;
                auto search_pos = exponential_search_upper_bound(childpieces,predicted_pos, key) -1;
                // 判断 找到的inner node 是否满足条件
                auto curit = (*(itlevel))->innerpieces.begin() + search_pos;
                FILMinsert<key_type,value_type>::Innerpiece* b = *curit;
//                if (b->startkey > key )
//                    cout << "i need You, my lovely Lord, b->startkey > key " << endl;

                for (; --itlevel >= innerlevels.begin();){
                    predicted_pos = (*b)(key);
                    childpieces = (*(itlevel))->innerpieces;
                    search_pos = exponential_search_upper_bound(childpieces,predicted_pos, key) -1;
                    // 判断 找到的inner node 是否满足条件
                    curit = (*(itlevel))->innerpieces.begin() + search_pos;
                }

                // find the leaf level
//                auto curleaf = leaflevel.leafpieces.begin();

                b = *curit;
                predicted_pos = (*b)(key);

                search_pos = exponential_search_upper_bound(leaflevel.leafpieces,predicted_pos, key) -1;
                auto curleaf = leaflevel.leafpieces.begin() + search_pos;

                /*
                Leafpiece* a = *curleaf;
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

                if (key != a->slotkey[resslot])
                {
                    // find in sort_list
                    sortPiece* buffer = (sortPiece*) a->buffer_piece;
                    auto sortpos = (*buffer)(key);
                    if (buffer->slotkey[sortpos] == key){
                        result_find index_res = result_find(false,buffer->slotflag[sortpos],buffer,sortpos);   // sort_list data
                        return index_res;
                    }
                    else{
                        cout << " i need You, my Lord, please help me! not find the key in sort piece" << endl;
                    }
                }
                result_find index_res = result_find(true,a->slotflag[resslot],a,resslot);   // regular data
                return index_res;
                */

                Leafpiece *a = *curleaf;
                // 确定 该leaf 的哪一个链 leaf piece
                while (key > a->endkey && a->buffer_piece) {
                    a = a->buffer_piece;
                }
                if (a->buffer_piece == NULL) {
                    sortPiece *buffer = (sortPiece *) a;
                    auto sortpos = (*buffer)(key);
                    if (buffer->slotkey[sortpos] == key) {
                        result_find index_res = result_find(false, buffer->slotflag[sortpos], buffer,
                                                            sortpos);   // sort_list data
                        return index_res;
                    } else {
                        cout << " i need You, my Lord, please help me! not find the key " << key << endl;
                    }
                }
                else {
                    predicted_pos = (*a)(key);

                    int slotlo = PGM_SUB_EPS2(predicted_pos, ErrorRecursive + 1, (a->slotkey.size() - 1));
                    int slothi = PGM_ADD_EPS(predicted_pos, ErrorRecursive, (a->slotkey.size() - 1));
                    int resslot = slotlo;
                    for (; resslot < slothi; resslot++) {
                        if (key <= a->slotkey[resslot]) {
                            break;
                        }
                    }

                    if (key != a->slotkey[resslot]) {
                        // find in sort_list
                        while (a->buffer_piece->buffer_piece != NULL) {
                            a = a->buffer_piece;
                        }
                        sortPiece *buffer = (sortPiece *) a->buffer_piece;
                        auto sortpos = (*buffer)(key);
                        if (buffer->slotkey[sortpos] == key) {
                            result_find index_res = result_find(false, buffer->slotflag[sortpos], buffer,
                                                                sortpos);   // sort_list data
                            return index_res;
                        } else {
                            cout << " i need You, my Lord, please help me! not find the key " << key << endl;
                        }
                    }

                    result_find index_res = result_find(true, a->slotflag[resslot], a, resslot);   // regular data
                    return index_res;
                }

            }
            else{

                // 只有一层leaf
                static constexpr size_t linear_search_threshold = 8 * 64 / sizeof(Innerpiece);   //sizeof(Innerpiece)

                // leaf level
                auto curleaf = leaflevel.leafpieces.begin();
                auto predicted_pos = (*root)(key);

                auto leaflo = PGM_SUB_EPS2(predicted_pos, Error + 1,leaflevel.leafpieces.size()) + leaflevel.leafpieces.begin();
                auto leafhi = PGM_ADD_EPS(predicted_pos, Error,leaflevel.leafpieces.size()-1) + leaflevel.leafpieces.begin();
                auto h = PGM_ADD_EPS(predicted_pos, Error,leaflevel.leafpieces.size()-1);
//
//                Leafpiece*  y= *leafhi;
//                Leafpiece*  z= *leaflo;
                if (key >= (*leafhi)->startkey){curleaf = leafhi;}
                else{
                    if constexpr (4 <= linear_search_threshold ){   //ErrorRecursive
                        for (; (*leaflo)->endkey < key; leaflo++){
//                            Leafpiece*  x= *leaflo;
                            continue;
                        }

                        curleaf = leaflo;
                    }
                    else{ curleaf = std::lower_bound(leaflo, leafhi, key);   // upper_bound 返回在前闭后开区间查找key 的上届，返回大于 key 的第一个元素位置； lower_bound 返回第一个大于或等于 key 的元素的位置
                    }
                }


                Leafpiece *a = *curleaf;
                // 确定 该leaf 的哪一个链 leaf piece
                while (key > a->endkey && a->buffer_piece) {
                    a = a->buffer_piece;
                }
                if (a->buffer_piece == NULL) {
                    sortPiece *buffer = (sortPiece *) a;
                    auto sortpos = (*buffer)(key);
                    if (buffer->slotkey[sortpos] == key) {
                        result_find index_res = result_find(false, buffer->slotflag[sortpos], buffer,
                                                            sortpos);   // sort_list data
                        return index_res;
                    } else {
                        cout << " i need You, my Lord, please help me! not find the key " << key << endl;
                    }
                }
                else {
                    predicted_pos = (*a)(key);

                    int slotlo = PGM_SUB_EPS2(predicted_pos, ErrorRecursive + 1, (a->slotkey.size() - 1));
                    int slothi = PGM_ADD_EPS(predicted_pos, ErrorRecursive, (a->slotkey.size() - 1));
                    int resslot = slotlo;
                    for (; resslot < slothi; resslot++) {
                        if (key <= a->slotkey[resslot]) {
                            break;
                        }
                    }

                    if (key != a->slotkey[resslot]) {
                        // find in sort_list
                        while (a->buffer_piece->buffer_piece != NULL) {
                            a = a->buffer_piece;
                        }
                        sortPiece *buffer = (sortPiece *) a->buffer_piece;
                        auto sortpos = (*buffer)(key);
                        if (buffer->slotkey[sortpos] == key) {
                            result_find index_res = result_find(false, buffer->slotflag[sortpos], buffer,
                                                                sortpos);   // sort_list data
                            return index_res;
                        } else {
                            cout << " i need You, my Lord, please help me! not find the key " << key << endl;
                        }
                    }

                    result_find index_res = result_find(true, a->slotflag[resslot], a, resslot);   // regular data
                    return index_res;
                }
            }
        }


        htaprange_res_find search_range(key_type firstkey, key_type lastkey,access_stats *query_stats) {
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
            auto leaf_res = index_res.findleafs[0]->find_upper_in_leaf(lastkey, Error);
            while (leaf_res.first == false) {
                Leafpiece *a = *curleaf;
                rangeresult.findleafs.push_back((*curleaf));
                rangeresult.leafslots.push_back(index_res.leafslot++);
                ++curleaf;
                leaf_res = (*curleaf)->find_upper_in_leaf(lastkey, Error);
            }
            rangeresult.findleafs.push_back((*curleaf));
            rangeresult.leafslots.push_back(index_res.leafslot++);
            rangeresult.slotleaf1 = leaf_res.second;
            int leaf_num= rangeresult.findleafs.size();
            if (leaf_num>1){
                sortPiece * buffer = (sortPiece *) (index_res.findleafs[1]);
                if (buffer->intercept >0)
                    rangeresult.findbuffers.emplace_back(buffer);

                // 检查新的 leaf 的buffer
                for (int i = 1; i<leaf_num;i++){
                    buffer = (sortPiece *) (rangeresult.findleafs[i]->buffer_piece);
                    while (buffer->buffer_piece != NULL){
                        Leafpiece * curleaf = (Leafpiece *)(buffer);
                        buffer = (sortPiece *) (curleaf->buffer_piece);
                    }
                    if (buffer->intercept > 0){   // = 0 说明buffer 中没有数据   446181678460  446345317732  446346181696
                        auto buffer_res = buffer->find_upper_in_leaf(lastkey, Error);
                        rangeresult.findbuffers.emplace_back(buffer);
                        rangeresult.slotbuffer1 = buffer_res.second;
                    }
                    else if(rangeresult.findbuffers.size()){
                        rangeresult.slotbuffer1 = rangeresult.findbuffers.back()->slotkey.size();
                    }

                }
            }
            else{
                sortPiece *buffer = (sortPiece *) (index_res.findleafs[1]);
                if (buffer->intercept > 0){   // = 0 说明buffer 中没有数据
                    auto buffer_res = buffer->find_upper_in_leaf(lastkey, Error);
                    rangeresult.findbuffers.emplace_back(buffer);
                    rangeresult.slotbuffer1 = buffer_res.second;
                }
            }


//            cout << "i turst in You, my Lord" << endl;
            return rangeresult;

        }


        template<class lru_type>
        forceinline bool insert_one(const key_type key, key_type* payload,lru_type interchain)  {
            // 判断一下，如果 小于 last leaf 的end key, 则插入到sortpiece 中
            if (key < this->m_tailleaf->endkey){
                // 找到属于哪个leaf piece
                // 插入到 该leaf piece 的buffer piece 中，
                // 判断插入到leaf piece 中，是否超出阈值
                auto search_res = this->search_one_leaf(key);
                Leafpiece* a = search_res.findleaf;
                int link_num = 0;
                while (a->buffer_piece->buffer_piece != NULL){
                    a = a->buffer_piece;
                    link_num ++;
                }
//                if (a->startkey == 446181678460 || a->startkey == 446345317732){
//                    cout << "O Lord, teach me where and how to seek You" <<endl;
//                }

                sortPiece* bufferpiece = (sortPiece*) (a->buffer_piece);
                auto flags = bufferpiece->insert(key,payload);

                if (flags.first)
                {  // findleaf rebuild
                    if (search_res.findleaf->startkey == 9615754926)
                        cout << "Jesus, i need You! please come" << endl;
                    key_type bufferkey = search_res.findleaf->learnedMergeSort(ErrorRecursive, interchain);
                    interchain->remove(bufferkey);  //
                    // retrain model, get the new slope and intercept
                    int newmodels = internal::make_segmentation( this,search_res.findleaf,Error,search_res.findleaf->slotkey,interchain,search_res.leafslot) - link_num;
                    this->newmodels += newmodels;
                    this->retrain ++;
//                    if (flags.second){
//                        cout << "i need You, my Lord, flags.second = true" <<endl;
//                    }
                }
                else{
//                    interchain->put((*search_res.findleaf).slotkey[0],search_res.findleaf);   // 如果要加入，需要判断， findleaf 的introchain.size 是否大于 1
                    if (flags.second){
                        interchain->remove( bufferpiece->slotkey[1]);
//                        cout << "i need You, my Lord" <<endl;
                    }
                    auto lrukey = bufferpiece->slotkey[0];
                    interchain->put(lrukey,bufferpiece);
                }

            }
            else{
                this->append_one(key, payload,Error, interchain);
            }
        }

        template<class lru_type>
        forceinline bool update_one(const key_type key, key_type* payload,lru_type interchain)  {
//  找到要 update 的key, 改变其payload。 原始payload 为 000000000， updated payload 为 1111111
            if (key < this->m_tailleaf->endkey){
                // 找到属于哪个leaf piece
                // 插入到 该leaf piece 的buffer piece 中，
                // 判断插入到leaf piece 中，是否超出阈值
                auto search_res = this->search_one_leaf(key);
                Leafpiece* a = search_res.findleaf;
                int link_num = 0;
                while (a->buffer_piece->buffer_piece != NULL){
                    a = a->buffer_piece;
                    link_num ++;
                }

                sortPiece* bufferpiece = (sortPiece*) (a->buffer_piece);
                auto flags = bufferpiece->insert(key,payload);


                if (flags.first)
                {  // findleaf rebuild
//                    if (search_res.findleaf->startkey == 60715691)
//                        cout << "Jesus, i need You! please come" << endl;
                    key_type bufferkey = search_res.findleaf->learnedMergeSort(ErrorRecursive, interchain);
                    interchain->remove(bufferkey);  //
                    // retrain model, get the new slope and intercept
                    int newmodels = internal::make_segmentation( this,search_res.findleaf,Error,search_res.findleaf->slotkey,interchain,search_res.leafslot) - link_num;
                    this->newmodels += newmodels;
                    this->retrain ++;
//                    if (flags.second){
//                        cout << "i need You, my Lord, flags.second = true" <<endl;
//                    }
                }
                else{
                    interchain->put((*search_res.findleaf).slotkey[0],search_res.findleaf);
                    if (flags.second){
                        interchain->remove( bufferpiece->slotkey[1]);
//                        cout << "i need You, my Lord" <<endl;
                    }
                    auto lrukey = bufferpiece->slotkey[0];
                    interchain->put(lrukey,bufferpiece);
                }

            }
            else{
                this->append_one(key, payload,Error, interchain);
            }
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


        const map<string, double> runtimeshow_verify()   // with the structure of reduced tracking usage
        {
            std::map<std::string,double> treeinfo;

            double innodesize = 24;  // startkey, double(slope,intercept)
            double leafnodesize = 32 ; // startkey, endkey, double(slope,intercept),  slotkey(算入了data usage 中，因为在introchain 中 只存储了value)
            double hashlrusize = 48;
            double no_hashsize = 20;  //prev,next, slot-short int
            double dataV = (valuesize+1)*8;
            double datausage = double(inkeynum*(dataV))/1048576;
            double addusage = double((exkeynum+inkeynum)*(1+8)+exkeynum*sizeof(key_type) )/1048576;  //
            // addusage: (exkeynum+inkeynum)*(1+8) —— bitmap + pointer;  exkeynum*8 —— pageid(ungined int),offset(unsigned int) ; exkeynum*sizeof(key_type) ——exkey;
            double indexusge = double((vstats.leaves*2+newmodels)*leafnodesize + vstats.innernum*innodesize + 8 )/1048576;
            // leaf nodes (leaf nodes * 2 due to the sorted buffer pieces, each leaf has a buffer piece,
            // newmodels indicates the number of leaf models linked in leaf level), inner nodes, root
            double lruusage = double(no_hashsize*inkeynum + hashlrusize * vstats.leaves)/1048576;
            treeinfo.insert(pair<std::string,int>("leaves",vstats.leaves));
            treeinfo.insert(pair<std::string,int>("inners",vstats.innernum));
            treeinfo.insert(pair<std::string,int>("levels",vstats.numlevel));
            treeinfo.insert(pair<std::string,double>("datausage",datausage));
            treeinfo.insert(pair<std::string,double>("addusage",addusage));  //key,flag,pageID,offset
            treeinfo.insert(pair<std::string,double>("indexusage",indexusge));
            treeinfo.insert(pair<std::string,double>("lruusage",lruusage));
            return treeinfo;
        }


        const map<string, double> show_verify()
        {
            std::map<std::string,double> treeinfo;
            this->verify_film(&vstats);
            double innodesize = sizeof(key_type) + 8+8;  // startkey, double(slope,intercept)
            double leafnodesize = 2 * sizeof(key_type) + 8 + 8 ; // startkey, endkey, double(slope,intercept),  slotkey(算入了data usage 中，因为在introchain 中 只存储了value)
            double hashlrusize = sizeof(key_type) + 8+ 16 + 16;
            double no_hashsize = 16+4;  //prev,next, slot-short int
            double dataV = (valuesize+1)*8;
            double datausage = double(inkeynum*(dataV))/1024/1024;
            double addusage = double((exkeynum+inkeynum)*(1+8)+exkeynum*sizeof(key_type) )/double(1048576);  // exkey, flag, pageID,offset   (pageid-int-8 offset short int), flag for data in memory
            // addusage: (exkeynum+inkeynum)*(1+8) —— bitmap + pointer;  exkeynum*8 —— pageid(ungined int),offset(unsigned int); exkeynum*sizeof(key_type) ——exkey;
            double indexusge = (vstats.leaves*2*leafnodesize+vstats.innernum*innodesize + 8 )/double(1048576);  //  leaf nodes, inner nodes, root
            double lruusage = double(no_hashsize*inkeynum + hashlrusize * vstats.leaves)/1024/1024;

            treeinfo.insert(pair<std::string,int>("leaves",vstats.leaves));
            treeinfo.insert(pair<std::string,int>("inners",vstats.innernum));
            treeinfo.insert(pair<std::string,int>("levels",vstats.numlevel));
            treeinfo.insert(pair<std::string,double>("datausage",datausage));
            treeinfo.insert(pair<std::string,double>("addusage",addusage));  //key,flag,pageID,offset
            treeinfo.insert(pair<std::string,double>("indexusage",indexusge));
            treeinfo.insert(pair<std::string,double>("lruusage",lruusage));
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

    template<typename key_type,typename value_type>
    struct FILMinsert<key_type,value_type>::Leafpiece{

        key_type startkey;
        key_type endkey;
        vector<key_type> slotkey;

        vector< bool > slotflag;   //adalru::Node<key_type, std::vector<key_type> >*
        vector< void*  > slotdata;   //adalru::Node<key_type, std::vector<key_type> >*
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


        inline int binary_search_upper_bound(vector<key_type> keys,int l, int r, const key_type key) const {
            while (l < r) {
                int mid = l + (r - l) / 2;
                auto psize = keys.size();
                auto piece = keys[mid];// by chao
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
                    key_type mkey = (*(lo+m-bound)) ;// chao
//        key_greater(mkey, key);
                    auto mkeydiff = mkey - key;  //chao
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
        key_type learnedMergeSort(int Error,lru_type interchain){

            auto curleaf = this;
            while(curleaf->buffer_piece->buffer_piece){
                curleaf = curleaf->buffer_piece;
//                auto s = this->slotkey.size();  // vect2.insert(vect2.begin(), vect1.begin(), vect1.end());
                copy( curleaf->slotkey.begin(), curleaf->slotkey.end(), back_inserter(this->slotkey));
                copy(curleaf->slotflag.begin(), curleaf->slotflag.end(), back_inserter(this->slotflag));
                copy(curleaf->slotdata.begin(), curleaf->slotdata.end(), back_inserter(this->slotdata));
                interchain->remove(curleaf->startkey);
            }
            sortPiece * bufferpiece = (sortPiece*) curleaf->buffer_piece;
            vector<key_type> slotkey;
            vector<bool> slotflag;
            vector<void*> slotdata;
            unsigned long slot = 0;
            // the first key insert use learned sort
            auto key = bufferpiece->slotkey[0];

            //            if (key == 31835864 || key == 31836190)
//                cout << "i need You, my lovely Lord, rebuild 30359458" << endl;
            unsigned long learnedPos = (*this)(key);
            // 根据learned pos 找到大于等于 buffer key 的第一个元素
            int slotlo = PGM_SUB_EPS2(learnedPos, Error + 1,(this->slotkey.size()-1)) ;  //,(--itlevel)->size()

            int resslot = slotlo;   // 这里不能再用 linear search 应该用 exponential search
            for (; resslot<this->slotkey.size();resslot++){
                if (key <= this->slotkey[resslot])
                {
//                    auto ku = this->slotkey[resslot];
//                    auto ku0 = this->slotkey[resslot-1];
                    break;
                }
            }

            resslot = exponential_search_upper_bound(this->slotkey,learnedPos, key);
            slotkey.insert(slotkey.end(),this->slotkey.begin()+slot,this->slotkey.begin()+resslot);
            slotkey.emplace_back(key);
            slotdata.insert(slotdata.end(),this->slotdata.begin()+slot,this->slotdata.begin()+resslot);
            slotdata.emplace_back(bufferpiece->slotdata[0]);
            slotflag.insert(slotflag.end(),this->slotflag.begin()+slot,this->slotflag.begin()+resslot);
            slotflag.emplace_back(bufferpiece->slotflag[0]);
            slot = resslot;

            for (int i = 1; i< bufferpiece->slope;i ++){
//                auto it = bufferpiece->slotkey.begin(); it !=bufferpiece->slotkey.end(); it++
                auto key = bufferpiece->slotkey[i];
//                if (key == 31835864 || key == 31836190 )
//                {
//                    auto kkkkkkkk = this->slotkey.size();
//                    cout << "i need You, my lovely Lord, rebuild 30359458" << endl;
//                }

// 依次比较 buffer 和 this->slotkey 中的元素，哪个小就插入到slotkey 中
                if (slot < this->slotkey.size()){
                    if (key < this->slotkey[slot]){
                        slotkey.emplace_back(key);
                        slotdata.emplace_back(bufferpiece->slotdata[i]);
                        slotflag.emplace_back(bufferpiece->slotflag[i]);
                    } else{
                        i --;
                        slotkey.emplace_back(this->slotkey[slot]);
                        slotdata.emplace_back(this->slotdata[slot]);
                        slotflag.emplace_back(this->slotflag[slot++]);
                    }
                }
                else{
                    slotkey.insert(slotkey.end(),bufferpiece->slotkey.begin()+i,bufferpiece->slotkey.begin()+bufferpiece->slotkey.size());
                    slotdata.insert(slotdata.end(),bufferpiece->slotdata.begin()+i,bufferpiece->slotdata.begin()+bufferpiece->slotdata.size());
                    slotflag.insert(slotflag.end(),bufferpiece->slotflag.begin()+i,bufferpiece->slotflag.begin()+bufferpiece->slotflag.size());
                    break;
                }
//                this->slotkey.insert(this->slotkey.begin()+resslot,key);
//                this->slotdata.insert(this->slotdata.begin()+resslot,bufferpiece->slotdata[i]);
//                this->slotflag.insert(this->slotflag.begin()+resslot,bufferpiece->slotflag[i]);
                // i 为元素对应的下标
            }
            auto nowkey = this->slotkey[slot];
            slotkey.insert(slotkey.end(),this->slotkey.begin()+slot,this->slotkey.end());
            slotdata.insert(slotdata.end(),this->slotdata.begin()+slot,this->slotdata.end());
            slotflag.insert(slotflag.end(),this->slotflag.begin()+slot,this->slotflag.end());

            auto kkkkk = slotkey.back();
            auto kkkkk1 = this->slotkey.back();

            this->slotkey = slotkey;
            this->slotdata = slotdata;
            this->slotflag = slotflag;

            bufferpiece->slope = 0;
            bufferpiece->intercept = 0;
//            cout << "i need You, my Lord" << endl;
//            要将bufferpiece 从interchain 也就是global LRU中删除
            interchain->remove(bufferpiece->slotkey[0]);  // 要将bufferpiece 从interchain 也就是global LRU中删除
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
            pair<bool,unsigned long> is_in(false,0);
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
                        auto upperkey = this->slotkey[slotlo];
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

    template<typename key_type,typename value_type>
    struct FILMinsert<key_type,value_type>::Innerpiece{

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

    template<typename key_type,typename value_type>
    struct FILMinsert<key_type,value_type>::Innerlevel{
        std::vector<FILMinsert<key_type,value_type>::Innerpiece* > innerpieces;
//        vector<internal::insertPWLF<key_type, int>*> opt;
        internal::insertPWLF<key_type, int>* opt;
        unsigned int pos = 0;
        unsigned int nextpos = 0;
        FILMinsert<key_type,value_type>::Innerpiece* i_innerpiece = new FILMinsert<key_type,value_type>::Innerpiece();  // 该innerlevel 的最后一个inner piece
    };

    template<typename key_type,typename value_type>
    struct FILMinsert<key_type,value_type>::Leaflevel{
        std::vector<FILMinsert<key_type,value_type>::Leafpiece* > leafpieces;
//        vector<internal::insertPWLF<key_type, int>*> opt;
        internal::insertPWLF<key_type, int>* opt;
        unsigned int pos = 0;
    };

#pragma pack(pop)

    typedef filminsert::FILMinsert< key_type, key_type* > filmadatype;

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

                    if (evictbuffer->intrachain.size <= 1){
                        filmmem->lru->remove(evictleaf->key);
                    }
                    auto evictslotV = evictbuffer->intrachain.poptail();

                    // the tail of the accessed leaf
                    gettimeofday(&lt2, NULL);
                    ltimeuse += (lt2.tv_sec - lt1.tv_sec) + (double) (lt2.tv_usec - lt1.tv_usec) / 1000000.0;
                    // reduced_evicttable
                    auto writeevict = filmmem->reduced_evictkeytoinpage(evictbuffer->slotkey[evictslotV->key], evictslotV->value,
                                                                        diskpage);
//                    auto writeevict = filmmem->evictkeytoinpage(evictbuffer->slotkey[evictslotV->key], evictslotV->value,
//                                                                diskpage);

                    evictbuffer->slotflag[evictslotV->key] = false;
                    evictbuffer->slotdata[evictslotV->key] = (void*)writeevict;
                    delete evictslotV;
                    evictslotV = NULL;
                }
                else{

                    auto evictslotV = evictleaf->value->intrachain.poptail();   // the tail of the accessed leaf

                    gettimeofday(&lt2, NULL);
                    ltimeuse += (lt2.tv_sec - lt1.tv_sec) + (double) (lt2.tv_usec - lt1.tv_usec) / 1000000.0;
                    auto writeevict = filmmem->reduced_evictkeytoinpage(evictleaf->value->slotkey[evictslotV->key], evictslotV->value,
                                                                        diskpage);
//                    filmmem->evictkeytoinpage(evictleaf->value->slotkey[evictslotV->key], evictslotV->value,
//                                              diskpage, writeevict);
                    evictleaf->value->slotflag[evictslotV->key] = false;
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
//                    if (evictbuffer->intrachain.size == 0){
//                        cout << "Jesus, please come~" <<endl;
//                    }
                    if (evictbuffer->intrachain.size <= 1){
                        filmmem->lru->remove(evictleaf->key);
                    }
                    auto evictslotV = evictbuffer->intrachain.poptail();

                    // the tail of the accessed leaf
                    gettimeofday(&lt2, NULL);
                    ltimeuse += (lt2.tv_sec - lt1.tv_sec) + (double) (lt2.tv_usec - lt1.tv_usec) / 1000000.0;
                    auto writeevict = filmmem->reduced_evictkeytoinpage(evictbuffer->slotkey[evictslotV->key], evictslotV->value,
                                                                        diskpage);
//                    auto writeevict = filmmem->evictkeytoinpage(evictbuffer->slotkey[evictslotV->key], evictslotV->value,
//                                                                diskpage);
                    evictbuffer->slotflag[evictslotV->key] = false;
                    evictbuffer->slotdata[evictslotV->key] = (void*)writeevict;
                    delete evictslotV;
                    evictslotV = NULL;
                }
                else{


                    auto evictslotV = evictleaf->value->intrachain.poptail();
                    // the tail of the accessed leaf
                    gettimeofday(&lt2, NULL);

                    auto writeevict = filmmem->reduced_evictkeytoinpage(evictleaf->value->slotkey[evictslotV->key], evictslotV->value,
                                                                        diskpage);
//                    auto writeevict = filmmem->evictkeytoinpage(evictleaf->value->slotkey[evictslotV->key], evictslotV->value,
//                                                                diskpage);

                    evictleaf->value->slotflag[evictslotV->key] = false;
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
        double ltimeuse = 0.0;

        if (filmmem->reused_index_pageid.size() != 0){
//            gettimeofday(&wdt1, NULL);
            auto wstart_time = std::chrono::high_resolution_clock::now();

            auto reused_pagenum = filmmem->reused_index_pageid.size();   // 判断是否有 reused page
            for (unsigned int i = 0; i< reused_pagenum; i++){
                // 读取reused page，并且一次写满该page
                unsigned int reused_pageid = filmmem->reused_pageid[i];
                unsigned long reused_index_pageid = filmmem->reused_index_pageid[i];
                // 为该page 创建一个inpage in memory
                auto reused_page = filmmem->createinmempage(reused_pageid,diskpage->pagesize,diskpage->recordnum);
                for (int k = 0; k < diskpage->recordnum; k ++){
                    // evict data, get the tail node, remove the tail of intrachain
//                    gettimeofday(&lt1, NULL);
                    auto lstart_time = std::chrono::high_resolution_clock::now();
                    auto evictleaf = filmmem->lru->get_tail1();
//                    if (evictleaf->key == 309575466846)
//                        cout << "Jesus, i need You" <<endl;
                    if (evictleaf->value->buffer_piece == NULL){
                        filmadatype::sortPiece * evictbuffer = (filmadatype::sortPiece*) (evictleaf->value);
                        if (evictbuffer->intrachain.size <= 1){
                            filmmem->lru->remove(evictleaf->key);
                        }
                        auto evictslotV = evictbuffer->intrachain.poptail(); // the tail of the accessed leaf

//                        gettimeofday(&lt2, NULL);
//                        ltimeuse += (lt2.tv_sec - lt1.tv_sec) + (double) (lt2.tv_usec - lt1.tv_usec) / 1000000.0;
                        ltimeuse += std::chrono::duration_cast<std::chrono::nanoseconds>(
                                std::chrono::high_resolution_clock::now() - lstart_time)
                                .count();
                        // reduced_evicttable
                        auto writeevict = filmmem->reduced_evictkeytoinpage(evictbuffer->slotkey[evictslotV->key], evictslotV->value,
                                                                            diskpage,reused_page, reused_index_pageid);

                        evictbuffer->slotflag[evictslotV->key] = false;
                        evictbuffer->slotdata[evictslotV->key] = (void*)writeevict;
                        delete evictslotV;
                        evictslotV = NULL;
                    }
                    else{
                        auto evictslotV = evictleaf->value->intrachain.poptail();   // the tail of the accessed leaf
                        ltimeuse += std::chrono::duration_cast<std::chrono::nanoseconds>(
                                std::chrono::high_resolution_clock::now() - lstart_time)
                                .count();   // ns
//                        ltimeuse += (lt2.tv_sec - lt1.tv_sec) + (double) (lt2.tv_usec - lt1.tv_usec) / 1000000.0;

                        auto writeevict = filmmem->reduced_evictkeytoinpage(evictleaf->value->slotkey[evictslotV->key], evictslotV->value,
                                                                            diskpage,reused_page, reused_index_pageid);

                        evictleaf->value->slotflag[evictslotV->key] = false;
                        evictleaf->value->slotdata[evictslotV->key] = (void*)writeevict;

                        delete evictslotV;
                        evictslotV = NULL;
                    }
                }

                if (!(reused_page->freespace))  //如果还能容纳一个record，那么就向该页写出，如果不能，那么就将当前的inpage写出disk，创建新的inpage
                {
                    filmmem->inmempages.emplace_back(reused_page);
                }
//                delete reused_page;
            }

            filmindex->inkeynum -= (diskpage->recordnum * reused_pagenum);
            filmindex->exkeynum += (diskpage->recordnum * reused_pagenum);
            vector<unsigned long >().swap(filmmem->reused_index_pageid);
            vector<pageid_type>().swap(filmmem->reused_pageid);

            // 需要单独将每个page 写出去磁盘，因为reused 的page 不一定连续
            int pagenum = filmmem->runtimeevictreusedpagestodisk(diskpage);
            r_stats->writenum += pagenum;
            wdtimeuse += std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::high_resolution_clock::now() - wstart_time)
                    .count();   // ns

//            gettimeofday(&wdt2, NULL);
//            wdtimeuse += (wdt2.tv_sec - wdt1.tv_sec) + (double) (wdt2.tv_usec - wdt1.tv_usec) / 1000000.0;
            wdtimeuse -= ltimeuse;
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
                    if (evictbuffer->intrachain.size <= 1){
                        filmmem->lru->remove(evictleaf->key);
                    }
                    auto evictslotV = evictbuffer->intrachain.poptail(); // the tail of the accessed leaf
//                    gettimeofday(&lt2, NULL);
//                    ltimeuse += (lt2.tv_sec - lt1.tv_sec) + (double) (lt2.tv_usec - lt1.tv_usec) / 1000000.0;
                    ltimeuse += std::chrono::duration_cast<std::chrono::nanoseconds>(
                            std::chrono::high_resolution_clock::now() - lstart_time)
                            .count();   // ns
                    auto writeevict = filmmem->reduced_evictkeytoinpage(evictbuffer->slotkey[evictslotV->key], evictslotV->value,
                                                                        diskpage);
                    evictbuffer->slotflag[evictslotV->key] = false;
                    evictbuffer->slotdata[evictslotV->key] = (void*)writeevict;
                    delete evictslotV;
                    evictslotV = NULL;
                }
                else{
//                    gettimeofday(&lt1, NULL);
                    auto evictslotV = evictleaf->value->intrachain.poptail();
                    // the tail of the accessed leaf
                    ltimeuse += std::chrono::duration_cast<std::chrono::nanoseconds>(
                            std::chrono::high_resolution_clock::now() - lstart_time)
                            .count();   // ns

//                    gettimeofday(&lt2, NULL);
//                    ltimeuse += (lt2.tv_sec - lt1.tv_sec) + (double) (lt2.tv_usec - lt1.tv_usec) / 1000000.0;
                    auto writeevict = filmmem->reduced_evictkeytoinpage(evictleaf->value->slotkey[evictslotV->key], evictslotV->value,
                                                                        diskpage);

                    evictleaf->value->slotflag[evictslotV->key] = false;
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
            wdtimeuse -= ltimeuse;
        }

        wdtimeuse = wdtimeuse/1e9;
        ltimeuse = ltimeuse/1e9;
        r_stats->wdisktime += wdtimeuse;
        r_stats->lrutime += ltimeuse;
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
            bool flag = curleaf->slotflag[k];
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

    void htap_range_is_in_or_out(int start, int end, filmadatype::Leafpiece* curleaf ,vector<pair<key_type,key_type*>> *totalres , pre_dict* dict,int sizevalue ,filmadamemory *filmmem,access_stats *r_stats){
        struct timeval lt1, lt2;
        double ltimeuse;
        for (int k = start;k < end;k++)
        {
            bool flag = curleaf->slotflag[k];
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
                unsigned long index_pageid = floor(int_reduced_evict/274)*274;
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

    void htap_range_is_in_or_out(int start, int end, filmadatype::sortPiece* curleaf ,vector<pair<key_type,key_type*>> *totalres , pre_dict* dict,int sizevalue ,filmadamemory *filmmem,access_stats *r_stats){
        struct timeval lt1, lt2;
        double ltimeuse;
        for (int k = start;k < end;k++)
        {
            bool flag = curleaf->slotflag[k];
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

                pair<pageid_type, pageOff_type>  reduced_writeevict (0,0);

                auto short_reduced_evict = (unsigned long) curleaf->slotdata[k];  //writeevict 指向的是 pageid and offset
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
                unsigned long index_pageid = floor(int_reduced_evict/274)*274;
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

    int htap_range_prepass(filmadamemory *filmmem,filmadatype::htaprange_res_find *index_result,vector<pair<key_type,key_type*>> *result,pre_dict* dict,filmadatype *filmindex,access_stats *r_stats){  // 将range query 由btree 得到的信息，进行怕热pass，遍历find leaf中的每一个leaf，如果slotdata.first = true, 则直接放入reslut 中，否则将其放入prepassdict 中

        int leafnum = index_result->findleafs.size();
        if (leafnum == 1){  //all the keys in range belonging to a certain leaf
//            cout << "Lord, this is chao" << endl;
            htap_range_is_in_or_out(index_result->slotleaf0,index_result->slotleaf1+1,index_result->findleafs[0],result,dict,filmindex->valuesize,filmmem,r_stats);
        }
        else {  // the first leaf and last leaf need to be deal with seperately
            // deal with firstleaf
//            cout << "Lord, this is chao" << endl;
            htap_range_is_in_or_out(index_result->slotleaf0,index_result->findleafs[0]->slotkey.size(),index_result->findleafs[0],result,dict,filmindex->valuesize,filmmem,r_stats);
            // deal with the middle leafs
            for (int i=1; i < leafnum-1; i++){
                htap_range_is_in_or_out(0,index_result->findleafs[i]->slotkey.size(),index_result->findleafs[i],result,dict,filmindex->valuesize,filmmem,r_stats);
            }
            //deal with lastleaf
//            cout << "Lord, this is chao" << endl;
            htap_range_is_in_or_out(0,index_result->slotleaf1+1,index_result->findleafs[leafnum-1],result,dict,filmindex->valuesize,filmmem,r_stats);
        }

        int buffernum = index_result->findbuffers.size();
        if (buffernum > 0){
            if (buffernum == 1){  //all the keys in range belonging to a certain leaf
//                cout << "Lord, this is chao" << endl;
                if (index_result->findbuffers[0]->intercept > 0)
                    htap_range_is_in_or_out(index_result->slotbuffer0,index_result->slotbuffer1,(filmadatype::sortPiece*)index_result->findbuffers[0],result,dict,filmindex->valuesize,filmmem,r_stats);
            }
            else {  // the first leaf and last leaf need to be deal with seperately
                // deal with firstleaf
//                cout << "Lord, this is chao" << endl;
                filmadatype::sortPiece* buffer = (filmadatype::sortPiece*) (index_result->findbuffers[0]);
                htap_range_is_in_or_out(index_result->slotbuffer0,buffer->slotkey.size(),(filmadatype::sortPiece*)index_result->findbuffers[0],result,dict,filmindex->valuesize,filmmem,r_stats);
                // deal with the middle leafs
                for (int i=1; i < buffernum-1; i++){
                    buffer = (filmadatype::sortPiece*) (index_result->findbuffers[i]);
                    htap_range_is_in_or_out(0,buffer->slotkey.size(),index_result->findbuffers[i],result,dict,filmindex->valuesize,filmmem,r_stats);
                }
                //deal with lastleaf
//                cout << "Lord, this is chao" << endl;
                htap_range_is_in_or_out(0,index_result->slotbuffer1,(filmadatype::sortPiece*)(index_result->findbuffers[buffernum-1]),result,dict,filmindex->valuesize,filmmem,r_stats);
            }
        }

        return 0;
    }


    // 一次至少读取1024k 的数据
    int range_read_from_disk(vector<pair<key_type,key_type*>> *result,pre_dict* dict,filmadamemory *filmmem,filmadatype *filmindex,filmadadisk *pagedisk, access_stats *r_stats) {
        // 打开磁盘文件
        struct timeval lt1, lt2, dt1, dt2, rdt1, rdt2;
        double ltimeuse = 0.0;  // 读磁盘的时间，写磁盘的时间
        double dtimeuse, rdtimeuse;
        int cross = dict->prepass.rbegin()->first - dict->prepass.begin()->first + 1;
//        cout<< cross << " ";
        r_stats->pagecross += dict->pagenum;
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
                unsigned long index_pageid = v.first*274;
                auto iiii = filmmem->reduced_evicttable[index_pageid];
//                if (pageid == 3027)
//                    cout << "Jesus, please teach me!" << endl;
                unsigned long int offset = diff * (pagedisk->recordsize * pagedisk->recordnum);// 把当前页的数据 单独出来
                for (int k = 0; k <
                                v.second.size(); k++) {   //v.second.size()  是要从当前页中读取多少个data， v.second[i] 是每个data 的相关信息是一个pair，v.second[i].first 是所属的leaf， v.second[i].second。first是slot in leaf，second 是 evictposi指向page id 和off
//                    auto slot1 = v.second[k].second.first;
                    auto slot = v.second[k].second.first;
                    pageOff_type off = v.second[k].second.second;   // offset in page
                    auto writeevict = v.second[k].second.second;


                    if (v.first == filmmem->inpage->pageid) {
                        //vector<key_type> inmemdata = filmmem->inpage->inmemdata;
                        key_type *inmemdata = filmmem->inpage->inmemdata;
                        res.first = inmemdata[off];
                        res.second = new key_type[filmindex->valuesize];
                        for (int i = 0; i < filmindex->valuesize; i++) {
                            res.second[i] = inmemdata[off + 1 + i];
                        }

                    } else {
                        res.first = buf[offset + off];
                        res.second = new key_type[filmindex->valuesize];
                        for (int i = 0; i < filmindex->valuesize; i++) {
                            res.second[i] = buf[offset + off + 1 + i];
                        }
                    }
                    auto recordnum = off/pagedisk->recordsize;  // compute the key's rank in the page
                    filmmem->unset_bit(recordnum, index_pageid);
                    filmmem->reduced_evicttable[index_pageid+1] -= 1;
                    filmindex->inkeynum++;
                    filmindex->exkeynum--;
                    result->push_back(res);

                    if (v.second[k].first->buffer_piece != NULL){
                        filmadatype::Leafpiece *curpiece = v.second[k].first;
                        if (res.first != curpiece->slotkey[slot]) {
                            cout << "i need You, my Lord, Thank You for all the thing! " << endl;
                        }
                        // 将res 插入到curpiece 的intrachain，// 将被访问的record 的数据放入内存中，即其所属的leaf，并修改其slotdata 的flag 为 true
                        gettimeofday(&lt1, NULL);
                        curpiece->slotdata[slot] = curpiece->intrachain.put(slot, res.second);
                        curpiece->slotflag[slot] = true;

                        // 将curpiece 插入到interchain 中，即mem->lru中
                        filmmem->lru->put(curpiece->startkey, curpiece);// 更新 interchain
                        gettimeofday(&lt2, NULL);
                        ltimeuse += (lt2.tv_sec - lt1.tv_sec) + (double) (lt2.tv_usec - lt1.tv_usec) / 1000000.0;
                    }
                    else{
                        filmadatype::sortPiece *curpiece = (filmadatype::sortPiece *)v.second[k].first;
                        if (res.first != curpiece->slotkey[slot]) {
                            cout << "i need You, my Lord, Thank You for all the thing! " << endl;
                        }
                        // 将res 插入到curpiece 的intrachain，// 将被访问的record 的数据放入内存中，即其所属的leaf，并修改其slotdata 的flag 为 true
                        gettimeofday(&lt1, NULL);
                        curpiece->slotdata[slot] = curpiece->intrachain.put(slot, res.second);
                        curpiece->slotflag[slot] = true;

                        // 将curpiece 插入到interchain 中，即mem->lru中
                        filmmem->lru->put(curpiece->startkey, curpiece);// 更新 interchain
                        gettimeofday(&lt2, NULL);
                        ltimeuse += (lt2.tv_sec - lt1.tv_sec) + (double) (lt2.tv_usec - lt1.tv_usec) / 1000000.0;
                    }

                }
                // check whether to merge key
                if (filmmem->check_merge(index_pageid) == true){
                    if (pageid != filmmem->inpage->pageid){
                        auto numinpage = filmmem->reduced_evicttable[index_pageid+1]; // record the index_pageid, that the position in the reduced_evictedtable, which can be reused.
                        filmmem->reused_pageid.emplace_back(pageid);
                        filmmem->num_page_merge += 1;
                        filmmem->reused_index_pageid.emplace_back(index_pageid);
                        int merge_num= filmmem->page_merge(buf,pagedisk->recordnum, index_pageid,offset);

                        filmindex->exkeynum -= filmmem->reduced_evicttable[index_pageid+1];
                        filmindex->inkeynum += filmmem->reduced_evicttable[index_pageid+1];
                        filmmem->reduced_evicttable[index_pageid+1] = 0;
                    }
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
                unsigned long index_pageid = pageid*274;
                auto iiii = filmmem->reduced_evicttable[index_pageid];
//                if (pageid == 3027)
//                    cout << "Jesus, plese teach me!" << endl;
                if (pageid == filmmem->inpage->pageid) {
                    cout << " Jesus, Son of David, please have pity on me, read from in-memory page" << endl;
                    key_type *inmemdata = filmmem->inpage->inmemdata;

                    for (int vi = 0; vi <
                                     v.second.size(); vi++) {   //v.second.size()  是要从当前页中读取多少个data， v.second[i] 是每个data 的相关信息是一个pair，v.second[i].first 是所属的leaf， v.second[i].second。first是slot in leaf，second 是 evictposi指向page id 和off
                        auto slot = v.second[vi].second.first;   // v.first 是 page id，v.second 是要从该页中读进来的每个record 的信息，是一个pair，
                        // 该pair 的first 是所属的leaf，该pair 的second 是一个pair，
                        // 包括：first 是在该leaf 中的slot，second 是page 信息，包括page id 和在该page 中的offset
                        short int off = v.second[vi].second.second;  //offset in page
                        auto writeevict = v.second[vi].second.second;
                        filmadatype::Leafpiece *curpiece = v.second[vi].first;
                        res.first = inmemdata[off];
                        res.second = new key_type[filmindex->valuesize];
                        for (int ki = 0; ki < filmindex->valuesize; ki++) {
                            res.second[ki] = inmemdata[off + 1 + ki];
                        }// check page merge
                        auto recordnum = off/pagedisk->recordsize;  // compute the key's rank in the page
                        bool former_flag = filmmem->check_exists(recordnum,index_pageid);
                        if (former_flag == false){
//                auto numinpage = memoryfilm->pageinfo[index_pageid+1];
                            cout << "i need You, my Lord, the flag is wrong" << endl;
                        }
                        filmmem->unset_bit(recordnum, index_pageid);
                        filmmem->reduced_evicttable[index_pageid+1] -= 1;
                        if (former_flag == false){
                            cout << "i need You, my Lord, the flag is wrong" << endl;
                        }

                        filmindex->inkeynum++;
                        filmindex->exkeynum--;
                        result->push_back(res);

                        if (res.first != curpiece->slotkey[slot]) {
                            cout << "i need You, my Lord, Thank You for all the thing! " << endl;
                        }
                        // 将res 插入到curpiece 的intrachain，// 将被访问的record 的数据放入内存中，即其所属的leaf，并修改其slotdata 的flag 为 true
                        gettimeofday(&lt1, NULL);
                        curpiece->slotdata[slot] = curpiece->intrachain.put(slot, res.second);
                        curpiece->slotflag[slot] = true;
                        // 将curpiece 插入到interchain 中，即mem->lru中
                        filmmem->lru->put(curpiece->startkey, curpiece);// 更新 interchain
                        gettimeofday(&lt2, NULL);
                        ltimeuse += (lt2.tv_sec - lt1.tv_sec) + (double) (lt2.tv_usec - lt1.tv_usec) / 1000000.0;
//                      cout<< " i need You, my Lord!!! we all need You! "<< endl;
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
                        filmadatype::Leafpiece *curpiece = v.second[vi].first;
                        res.first = buf[off];
                        res.second = new key_type[filmindex->valuesize];
                        for (int ki = 0; ki < filmindex->valuesize; ki++) {
                            res.second[ki] = buf[off + 1 + ki];
                        }
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
                            curpiece->slotdata[slot] = curpiece->intrachain.put(slot, res.second);
                            curpiece->slotflag[slot] = true;
                            // 将curpiece 插入到interchain 中，即mem->lru中
                            filmmem->lru->put(curpiece->startkey, curpiece);// 更新 interchain
                            gettimeofday(&lt2, NULL);
                            ltimeuse += (lt2.tv_sec - lt1.tv_sec) + (double) (lt2.tv_usec - lt1.tv_usec) / 1000000.0;
//                cout<< " i need You, my Lord!!! we all need You! "<< endl;
                            if (res.first != curpiece->slotkey[slot]) {
                                cout << "i need You, my Lord, Thank You for all the thing! " << endl;
                            }
                        }
                        else{
                            filmadatype::sortPiece *curpiece = (filmadatype::sortPiece *)v.second[vi].first;
                            // 将res 插入到curpiece 的intrachain，// 将被访问的record 的数据放入内存中，即其所属的leaf，并修改其slotdata 的flag 为 true
                            gettimeofday(&lt1, NULL);
                            curpiece->slotdata[slot] = curpiece->intrachain.put(slot, res.second);
                            curpiece->slotflag[slot] = true;
                            // 将curpiece 插入到interchain 中，即mem->lru中
                            filmmem->lru->put(curpiece->startkey, curpiece);// 更新 interchain
                            gettimeofday(&lt2, NULL);
                            ltimeuse += (lt2.tv_sec - lt1.tv_sec) + (double) (lt2.tv_usec - lt1.tv_usec) / 1000000.0;
                            if (res.first != curpiece->slotkey[slot]) {
                                cout << "i need You, my Lord, Thank You for all the thing! " << endl;
                            }
//                cout<< " i need You, my Lord!!! we all need You! "<< endl;
                        }
                    }
                }

                if (filmmem->check_merge(index_pageid) == true){

                    if (pageid != filmmem->inpage->pageid) {
                        auto numinpage = filmmem->reduced_evicttable[index_pageid +1]; // record the index_pageid, that the position in the reduced_evictedtable, which can be reused.
                        filmmem->reused_pageid.emplace_back(pageid);
                        filmmem->num_page_merge += 1;
                        filmmem->reused_index_pageid.emplace_back(index_pageid);
                        int merge_num = filmmem->page_merge(buf, pagedisk->recordnum, index_pageid);

//                    if (numinpage != merge_num)
//                        cout<< "Jesus, i trust in You! numinpage != merge_num" << endl;
                    filmindex->exkeynum -= filmmem->reduced_evicttable[index_pageid+1];
                    filmindex->inkeynum += filmmem->reduced_evicttable[index_pageid+1];
                    filmmem->reduced_evicttable[index_pageid+1] = 0;
                    }
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
                                     unsigned long init_num_keys, double insert_frac,double out_of_order_frac, unsigned long batch_size,string workload,unsigned int  rebuild_threshold) {

        filmadatype filmada(0, errbnd, errbnd);
        filmadalrutype *interchain = new filmadalrutype(actual_numkey);
//        auto kkk = int (actual_numkey/600);
        filmada.leaflevel.leafpieces.reserve(int (actual_numkey/600));
        filmada.valuesize = recordSize - 1;

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
        filmstorage::filmmemory<key_type, key_type *, filmadatype, filmadalrutype, filmadadisk> memoryfilm(numkey,
                                                                                                           mem_threshold,
                                                                                                           &filmada,
                                                                                                           interchain);
        memoryfilm.reserveMem = reserveMem;
        memoryfilm.merge_threshold = merge_threshold;
        filmstorage::filmdisk<key_type> diskfilm(diskfile, pagesize, numrecord, recordSize, merge_threshold);

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

        key_type payload[filmada.valuesize]{};
        key_type update_payload[filmada.valuesize]{1};

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
        string performance_file = "/home/wamdm/chaohong/clionDir/updatefilm/result/totally_rebuild_reduced_adalru_film_performance.txt";
        savefile.open(performance_file, ios::app);
        savefile << "method " << "film_ada_lru " << "available_memory " << mem_threshold + reserveMem << " qtype "
                 << workload << " error " << errbnd << " pagesize " << (pagesize * 8 / 1024) ;
        savefile << "k recordsize "<< recordSize << " build_time " << buildtimeuse << " dataset " << dataset << " datasize " << datasize
                 << " keynum " << numkey << " merge_threshold " << merge_threshold <<" rebuild_threshold " << rebuild_threshold << " workload " <<workload <<" ";
        savefile << endl;
        savefile << flush;
        savefile.close();

        savefile.open(performance_file, ios::app);
        map<string, double>::iterator iter;
        savefile << "init_state " << "film ";
        savefile << "init_write_time " << diskfilm.initwtime / 1000000.0  << " ";
        savefile << "pages_init_write " << diskfilm.nextpageid << " ";
        savefile << "init_num_keys " << init_num_keys << " ";
        savefile << "method " << "film_htap_flatten ";
        savefile << "out_of_order_frac "<< out_of_order_frac << " ";
        for (iter = transflag.second.meminfo.begin(); iter != transflag.second.meminfo.end(); iter++)
            savefile << iter->first << " " << iter->second << " ";
        savefile << "\n";
        savefile << flush;
        savefile.close();

        // Run workload
        unsigned int i = init_num_keys;
        long long cumulative_inserts = 0, cumulative_updates = 0, retrain_updates = 0,
            cumulative_appends = 0, cumulative_ranges = 0;
        long long cumulative_lookups = 0;

        double cumulative_insert_time = 0;
        double comulative_rebuild_time = 0;
        double cumulative_lookup_time = 0;
        double cumulative_update_time = 0, cumulative_append_time = 0;
        double cumulative_range_time = 0;
        string lookup_distribution = workload;
        struct timeval workload_start_time, workload_run_time, insert_start_time, insert_end_time,
                lookup_start_time, lookup_end_time;
        double workload_elapsed_time = 0.0, lookup_elapsed_time = 0.0, insert_elapsed_time = 0.0,
            update_elapsed_time = 0.0, append_elapsed_time = 0.0,rebuild_elapsed_time = 0.0,
            range_elapsed_time = 0;
        double time_limit = 7200;
//        double querycumulative_writetime = 0.0;
        double cumulative_query_writetime = 0.0;
        double cumulative_insert_writetime = 0; // 记录插入process总共导致了多少 写磁盘的时间
        gettimeofday(&workload_start_time, NULL);
        unsigned int batch_no = 0;
        unsigned int total_retrain = 0;
        int total_newmodels = 0;


        access_stats query_stats;
        bool print_batch_stats = true;
        query_stats.qtype = "query";
        query_stats.zipffactor = zipf;
        query_stats.workload = workload;
        query_stats.insert_frac = insert_frac;
        query_stats.out_of_order_frac = out_of_order_frac ;

        pair<bool, double> runtime_transflag = memoryfilm.runtimejudgetrans(&query_stats);
        int periodV = memoryfilm.reserveMem * 1024 * 1024 / ((filmada.valuesize + 1) * sizeof(key_type) * 10);
        enum Operation {
            READ = 0, INSERT, DELETE, SCAN, UPDATE, DIFF
        };
        unsigned long int num_keys_after_batch = work.init_table_size;
        while (true) {
            batch_no++;
//            if (batch_no == 83 )
//                cout << "i need You, my Lord" << endl;
            // generate workload according to the work.parameters
            // generate the search key within the inserted key space

            if (i > 0) {
                work.generate_operations(work.keys, num_keys_after_batch);

                // step1: do lookup in htap;
                query_stats.computetimeuse = 0.0;
                double batch_query_writetime = 0.0;
                unsigned long num_actual_lookups = work.read_ops.size();
                if (num_actual_lookups) {
                    int qi = 0;
                    gettimeofday(&lookup_start_time, NULL);
                    for (auto &op_key: work.read_ops) {
//                         if (op_key == 291660467280 ) {
//                            cout << "my Lord, i need You! thank You, for ever and ever ~~" << endl;
//                        }
                        pair<key_type, key_type *> res = memoryfilm.get_key(op_key, &query_stats, diskfilm);
                        if (query_stats.disknum > 0 && qi++ % periodV == 0) {
                            runtime_transflag = memoryfilm.runtimejudgetrans(&query_stats);
                            while (runtime_transflag.first) {
                                batch_query_writetime += runtimeevictkeytopage2(&memoryfilm, runtime_transflag.second,
                                                                                &filmada, &diskfilm, &query_stats);
                                runtime_transflag = memoryfilm.runtimejudgetrans(&query_stats);
                            }
                        }
//                        cout << "this is read, thank You, my Lord" << endl;
                    }

                    // record the start time of lookup
                    gettimeofday(&lookup_end_time, NULL);
                    lookup_elapsed_time = (lookup_end_time.tv_sec - lookup_start_time.tv_sec) +
                                          (double) (lookup_end_time.tv_usec - lookup_start_time.tv_usec) / 1000000.0;
                    lookup_elapsed_time -= query_stats.computetimeuse;
                    cumulative_lookup_time += lookup_elapsed_time;
                    cumulative_lookups += num_actual_lookups;
                    cumulative_query_writetime += batch_query_writetime;
                    vector<key_type>().swap(work.read_ops);
                }


                // step2: Do inserts, with htap

                unsigned last_retrain = filmada.retrain;
                int last_newmodels = filmada.newmodels;
                double batch_insert_writetime = 0;   // sec
                unsigned long int num_actual_inserts = work.insert_ops.size();
                cout<< "batch " << batch_no << "num_actual_inserts: " << num_actual_inserts;
                int batch_newmodels = 0;
                unsigned batch_retrain = 0;
                if (num_actual_inserts){
                    gettimeofday(&insert_start_time, NULL);   //
                    auto start_time = std::chrono::high_resolution_clock::now();

                    filmada.insert_random(work.insert_ops, payload, interchain);
                    filmada.root = filmada.innerlevels.back()->innerpieces[0];
                    auto a = filmada.leaflevel.opt->get_segment(filmada.m_tailleaf->endkey);
                    filmada.m_tailleaf->update(a);

                    insert_elapsed_time =
                    std::chrono::duration_cast<std::chrono::nanoseconds>(
                            std::chrono::high_resolution_clock::now() - start_time)
                            .count();
//                    gettimeofday(&insert_end_time, NULL);
//                    insert_elapsed_time = (insert_end_time.tv_sec - insert_start_time.tv_sec) +
//                                          (double) (insert_end_time.tv_usec - insert_start_time.tv_usec) / 1000000.0;
//                    cout << "this is insert, thank You, my Lord" << endl;

                    // step3: judge whether to evict data, after doing inserts
                    transflag = memoryfilm.judgetransfer(&query_stats);
                    query_stats.computetimeuse = 0.0;
                    while (transflag.first) {
//                memoryfilm.runtimefilmtransfer(&diskfilm);
                        batch_insert_writetime += runtimeevictkeytopage2(&memoryfilm, transflag.second.totalusemem,
                                                                         &filmada, &diskfilm, &query_stats);
                        transflag = memoryfilm.judgetransfer(&query_stats);
                    }

                    vector<key_type>().swap(work.insert_ops);

                    batch_newmodels = filmada.newmodels - last_newmodels;
                    batch_retrain = filmada.retrain - last_retrain;
                    total_retrain += batch_retrain;
                    total_newmodels += batch_newmodels;

                    if (filmada.newmodels >= rebuild_threshold){
                        gettimeofday(&lookup_start_time, NULL);
                        auto start_time = std::chrono::high_resolution_clock::now();
                        // rebuild internal levels
                        filmada.internal_level_rebuild(errbnd,interchain);
                        filmada.newmodels = 0;
                        filmada.retrain = 0;
                        filmada.rebuild_inner += 1;

                        rebuild_elapsed_time = std::chrono::duration_cast<std::chrono::nanoseconds>(
                                std::chrono::high_resolution_clock::now() - start_time)
                                .count();
//                        gettimeofday(&lookup_end_time, NULL);
//                        rebuild_elapsed_time = (lookup_end_time.tv_sec - lookup_start_time.tv_sec) +
//                                               (double) (lookup_end_time.tv_usec - lookup_start_time.tv_usec) / 1000000.0;
                    }

                    // 将写出磁盘的时间算到插入时间里，更新插入的时间
                    insert_elapsed_time += rebuild_elapsed_time;
                    comulative_rebuild_time += rebuild_elapsed_time;
                    insert_elapsed_time += batch_insert_writetime;
//                insert_elapsed_time -= query_stats.computetimeuse;

                    cumulative_inserts += num_actual_inserts;
                    num_keys_after_batch += num_actual_inserts;
//                cumulative_appends += alen;
                    cumulative_insert_time += insert_elapsed_time;
                    cumulative_insert_writetime += batch_insert_writetime;
                }


                // step 4 do update, change the payload
                unsigned long int num_actual_updates = work.update_ops.size();
                if (num_actual_updates){
                    gettimeofday(&insert_start_time, NULL);   //
                    for (auto &op_key:work.update_ops) {
                        if (op_key == 291660467280 ) {
                            cout << "my Lord, i need You! thank You, for ever and ever ~~" << endl;
                        }
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
                unsigned long int num_actual_ranges = work.scan_ops.size();
                if (work.scan_ops.size()){
                    gettimeofday(&insert_start_time, NULL);   //
//                    for (auto &op_key:work.scan_ops) {
////                        pair<key_type, key_type *> res = memoryfilm.get_scan(op_key, &query_stats,diskfilm);
//                        cout << "this is scan, thank You, my Lord" << endl;
//                    }
                    for (unsigned int ri = 0; ri < num_actual_ranges; ri++) {
                        vector<pair<key_type, key_type *>> rres;
                        key_type firstkey = work.range_ops[ri][0];
                        key_type lastkey = work.range_ops[ri][1];
                        if (firstkey == 42265599683 || firstkey == 42264106012)   // 42266346608
                            cout << "Lord, please teach me to seek You" << endl;
                        auto index_res = memoryfilm.get_range(firstkey,lastkey, &query_stats,diskfilm);
                        pre_dict *dict = new pre_dict();
                        htap_range_prepass(&memoryfilm, &index_res, &rres, dict, &filmada, &query_stats);

                        // determine whether to access disk
                        if (dict->pagenum == 0) {  // the request data  数据都在内存
                            //          cout<< "Jesus, happy birthday!" << endl;
                            query_stats.memnum += 1;
                            // 遍历rres，释放数组
                        } else if (rres.empty()) { //数据都在磁盘
                            range_read_from_disk(&rres, dict, &memoryfilm, &filmada, &diskfilm,
                                                 &query_stats);// read data from disk according to 根据prepass 的信息，从磁盘中读数据

                            query_stats.disknum += 1;
                            query_stats.diskpagenum += dict->pagenum;
                            //cout << "i want in Your heart, my Lord!" << endl;
                        } else {   //the request data are in both memory and disk
                            query_stats.crossnum += 1;
                            query_stats.crosspagenum += dict->pagenum;
                            //                cout<< "Jesus, sister nana needs You!"<<endl;
                            range_read_from_disk(&rres, dict, &memoryfilm, &filmada, &diskfilm, &query_stats);
                        }
                        delete dict;
//                        cout << "this is range query, thank You, my Lord" << endl;
                    }
                    gettimeofday(&insert_end_time, NULL);
                    range_elapsed_time = (insert_end_time.tv_sec - insert_start_time.tv_sec) +
                                          (double) (insert_end_time.tv_usec - insert_start_time.tv_usec) / 1000000.0;
                    cumulative_range_time += range_elapsed_time;
                    cumulative_ranges += num_actual_ranges;
                   for (unsigned int ri = 0; ri < num_actual_ranges; ri++)
                        delete[] work.range_ops[ri];
                    delete[] work.range_ops;
                    work.range_ops = nullptr;
                    vector<key_type>().swap(work.scan_ops);
                }
                if (work.delete_ops.size()){
                    for (auto &op_key:work.read_ops) {
                        pair<key_type, key_type *> res = memoryfilm.get_key(op_key, &query_stats,diskfilm);
                        cout << "this is delete, thank You, my Lord" << endl;
                    }
                    num_keys_after_batch -= work.delete_ops.size();
                    vector<key_type>().swap(work.delete_ops);
                }


                // the total time
                gettimeofday(&workload_run_time, NULL);
                workload_elapsed_time = (workload_run_time.tv_sec - workload_start_time.tv_sec) +
                                        (double) (workload_run_time.tv_usec - workload_start_time.tv_usec) / 1000000.0;
//            if (workload_elapsed_time > time_limit && (filmada.inkeynum + filmada.exkeynum) > actual_numkey ) {
//                break;
//            }

                savefile.open(performance_file, ios::app);
                map<string, double>::iterator iter;
//                savefile << "_ ";
                savefile << "init_write_time " << diskfilm.initwtime / 1000000.0  << " ";
                savefile << "pages_init_write " << diskfilm.nextpageid << " ";
                savefile << "method " << "film_ada_lru ";
                savefile << "rebuild_threshold " << rebuild_threshold << " ";
                for (iter = transflag.second.meminfo.begin(); iter != transflag.second.meminfo.end(); iter++)
                    savefile << iter->first << " " << iter->second << " ";
//                savefile << "\n";
                savefile << flush;
                savefile.close();
                query_stats.print_stats();

                // step 5: print 输出
                if (print_batch_stats) {
                    int num_batch_operations = num_actual_lookups + num_actual_inserts +
                            num_actual_updates + num_actual_ranges;
                    double batch_time = lookup_elapsed_time + insert_elapsed_time +
                            update_elapsed_time + range_elapsed_time;
                    long long cumulative_operations = cumulative_lookups + cumulative_inserts +
                            cumulative_updates + cumulative_ranges;
                    double cumulative_time = cumulative_lookup_time + cumulative_insert_time +
                            cumulative_update_time + cumulative_range_time;

                    std::cout << "Batch " << batch_no
                              << ", cumulative_ops: " << cumulative_operations
                              << "\n\tbatch_throughput:\t"
                              << num_actual_lookups / lookup_elapsed_time
                              << " lookups/sec,\t"
                              << num_actual_inserts / insert_elapsed_time
                              << " inserts/sec,\t"
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
                              << cumulative_updates / cumulative_update_time
                              << " updates/sec,\t"
                            << cumulative_ranges / cumulative_range_time
                            << " ranges/sec,\t"
                              << cumulative_operations / cumulative_time << " ops/sec "
                              << cumulative_query_writetime << " cumulative_query_writetime "
                              << cumulative_insert_writetime << " cumulative_insert_writetime "
                              << rebuild_elapsed_time << " rebuild_elapsed_time "
                                << update_elapsed_time << " update_elapsed_time ";

                    std::cout << "\n\tcumulative_execution_1:\t"
                              << (filmada.inkeynum + filmada.exkeynum) * (filmada.valuesize + 1) * 8 /
                                 double(1024 * 1024)
                              << " G datasize, \t"
                              << filmada.inkeynum + filmada.exkeynum << " #cumulative_keys,\t"
                              << filmada.inkeynum << " #inmem_keys,\t"
                              << filmada.exkeynum << " #exmem_keys"
                              << "\n\tcumulative_execution_1:\t"
                              << cumulative_updates << " #cumulative_updates,\t" << cumulative_update_time
                              << " cumulative_update_time,\t"
                              << cumulative_inserts << " #cumulative_inserts,\t" << cumulative_insert_time
                              << " cumulative_insert_time,\t"
                              << cumulative_lookups << " #cumulative_lookups,\t" << cumulative_lookup_time
                              << " cumulative_lookup_time,\t"
                            << cumulative_ranges << " #cumulative_ranges,\t" << cumulative_range_time
                            << " cumulative_range_time,\t"
                             << total_retrain << " cumulative_retrain,\t "
                              <<  comulative_rebuild_time<< " comulative_rebuild_time,\t"
                              << (diskfilm.initwtime / 1000000.0) + query_stats.wdisktime << " cumulative_writetime "
                              << cumulative_insert_writetime << " cumulative_insert_writetime "
                              << cumulative_update_time << " cumulative_update_time ";

                    std::cout << "\n\tfilm_update_1:\t"
                              << rebuild_elapsed_time << " rebuild_elapsed_time "
                              << memoryfilm.num_page_merge << " num_page_merge "
                              << memoryfilm.merge_timeuse << " merge_timeuse "
                              << total_retrain << " times_of_leaves_retrain "
                              << total_newmodels << " updated_new_models "
                              << batch_retrain<< " batch_times_of_leaves_retrain "
                              << batch_newmodels << " batch_updated_new_models ";

                    std::cout << std::endl;

                    savefile.open(performance_file,ios::app);

                    savefile << batch_no << " Batch " <<
                             cumulative_operations << " cumulative_ops "
                             << "_" << " batch_throughput "
                             << num_actual_lookups / lookup_elapsed_time << " lookups/sec "
                             << num_actual_inserts / insert_elapsed_time << " inserts/sec "
                             << num_actual_updates / update_elapsed_time << " updates/sec "
                            << num_actual_ranges / range_elapsed_time << " ranges/sec "
                             << num_batch_operations / batch_time << " ops/sec "
                             << batch_query_writetime << " batch_query_writetime "
                             << batch_insert_writetime << " batch_insert_writetime "
                             << "_" << " cumulative_throughput "
                             << cumulative_lookups / cumulative_lookup_time << " cum_lookups/sec "
                             << cumulative_inserts / cumulative_insert_time << " cum_inserts/sec "
                             << cumulative_updates / cumulative_update_time << " cum_updates/sec "
                            << cumulative_ranges / cumulative_range_time << " cum_rangees/sec "
                             << cumulative_operations / cumulative_time << " cum_ops/sec "
                             << cumulative_query_writetime << " cumulative_query_writetime "
                             << cumulative_insert_writetime << " cumulative_insert_writetime "
                             << rebuild_elapsed_time << " rebuild_elapsed_time "
                             << filmada.rebuild_inner << " num_of_inner_rebuild "
                             << memoryfilm.num_page_merge << " num_page_merge "
                             << memoryfilm.merge_timeuse << " merge_timeuse ";


                    savefile << (filmada.inkeynum + filmada.exkeynum) * (filmada.valuesize + 1) * 8 / double(1024 * 1024)
                             << "G datasize "
                             << filmada.inkeynum + filmada.exkeynum << " #cumulative_keys "
                             << filmada.inkeynum << " #inmem_keys "
                             << filmada.exkeynum << " #exmem_keys "
                             << cumulative_lookups << " #cumulative_lookups " <<
                             cumulative_lookup_time  << " cumulative_lookup_time "
                            << cumulative_inserts << " #cumulative_inserts "
                            << cumulative_insert_time << " cumulative_insert_time "
                            << cumulative_updates << " #cumulative_updates "
                            << cumulative_update_time<< " cumulative_update_time "
                            << cumulative_ranges << " #cumulative_ranges "
                            << cumulative_range_time<< " cumulative_range_time "
                             << (diskfilm.initwtime / 1000000.0) + query_stats.wdisktime << " cumulative_writetime "
                             << memoryfilm.num_page_merge << " num_page_merge "
                             << memoryfilm.merge_timeuse << " merge_timeuse "
                             << total_retrain << " times_of_leaves_retrain "
                             << total_newmodels << " updated_new_models "
                             << batch_retrain<< " batch_times_of_leaves_retrain "
                             << batch_newmodels << " batch_updated_new_models ";

                    savefile << "\n";
                    savefile << flush;
                    savefile.close();

                }


                // Check for workload end conditions
                if (workload_elapsed_time > time_limit || (filmada.inkeynum + filmada.exkeynum) >= actual_numkey) {
                    break;
                }

//                if (num_actual_inserts < num_inserts_per_batch) {
//                    // End if we have inserted all keys in a workload with inserts
//                    break;
//                }
            }
            cout << "my lovely Lord, finished the interleave inserts and queries of batch " << batch_no << endl;

        }
        fs.open(diskfile,ios::in);
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
index_res= result_find(false,a->slotflag[resslot],a,resslot);   // not find, then, insert into bufferpiece,
index_res.leafslot = search_pos;
}
else{
index_res = result_find(true,a->slotflag[resslot],a,resslot);   // regular data
index_res.leafslot = search_pos;
}

return index_res;
 */