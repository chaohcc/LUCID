//
// Created by CCMa on 2022/1/28.
//

#ifndef FILMINSERT_PWLF_H
#define FILMINSERT_PWLF_H


#include <cmath>
#include <limits>
#include <vector>
#include <stdexcept>
#include <type_traits>
#include <climits>
#include <malloc.h>

#include <omp.h>
#include <iostream>
#include "film_plus.h"
#include "filmadastorage.h"

//#define forceinline inline __attribute__((__always_inline__))
using namespace std;

namespace filminsert::internal{
//    typedef filminsert::Leafpiece sort_piece_type;
//    std::conditional_t< 条件std::is_floating_point_v<T>, 条件为真 则 long double,条件为假 则 std::conditional_t<(sizeof(T) < 8), int64_t, __int128>>;
    template<typename T>    // std::is_floating_point_v<T>  检查T 是否为浮点类型，是则为true， else false
    using LargeSigned = typename std::conditional_t<std::is_floating_point_v<T>,
            long double,
            std::conditional_t<(sizeof(T) <8), int64_t, __int128>>;  // chaochao modify  "< " to <=

    template<typename X, typename Y>
    class insertPWLF{
        using SX = LargeSigned<X>;
        using SY = LargeSigned<Y>;
        struct Slope {
            SX dx{};
            SY dy{};

            bool operator<(const Slope &p) const { return dy * p.dx < dx * p.dy; }  // operator 运算符重载
            bool operator>(const Slope &p) const { return dy * p.dx > dx * p.dy; }
            bool operator==(const Slope &p) const { return dy * p.dx == dx * p.dy; }
            bool operator!=(const Slope &p) const { return dy * p.dx != dx * p.dy; }
            explicit operator long double() const { return dy / (long double) dx; }
        };

        struct Point {
            X x{};
            Y y{};

            Slope operator-(const Point &p) const { return {SX(x) - p.x, SY(y) - p.y}; }
        };

    public:
        const Y epsilon;
        std::vector<Point> lower;
        std::vector<Point> upper;
        X first_x = 0;
        X last_x = 0;
        unsigned long lower_start = 0;
        unsigned long  upper_start = 0;
        unsigned long points_in_hull = 0;
        Point rectangle[4];
        queue< Point > buffqueue;



        auto cross(const Point &O, const Point &A, const Point &B) const {
            auto OA = A - O;
            auto OB = B - O;

//        long int tmp = OA.dx * OB.dy - OA.dy * OB.dx;
//        std::cout<< tmp <<std::endl;
            return OA.dx * OB.dy - OA.dy * OB.dx;
        }

    public:

        class CanonicalSegment;

        explicit insertPWLF(Y epsilon) : epsilon(epsilon), lower(), upper() {
            if (epsilon < 0)
                throw std::invalid_argument("epsilon cannot be negative");

            upper.reserve(1u << 16);  // reserve, 告知容器应该准备保存多少元素，不改变容器中元素的数量，仅影响容器预先分配多少的内存空间
            lower.reserve(1u << 16);

        }
        explicit insertPWLF() : epsilon(), lower(), upper() {

            upper.reserve(1u << 16);  // reserve, 告知容器应该准备保存多少元素，不改变容器中元素的数量，仅影响容器预先分配多少的内存空间
            lower.reserve(1u << 16);
        }
        ~insertPWLF(){
            std::vector<Point>().swap(upper);
            std::vector<Point>().swap(lower);
//            delete [] rectangle;
            queue< Point >().swap(buffqueue);
            malloc_trim(0);
        }

        bool add_point(const X &x, const Y &y) {
//        if (points_in_hull > 0 && x <= last_x)
//            throw std::logic_error("Points must be increasing by x.");

            last_x = x;
            auto max_y = std::numeric_limits<Y>::max();
            auto min_y = std::numeric_limits<Y>::lowest();
            Point p1{x, y >= max_y - epsilon ? max_y : y + epsilon};
            Point p2{x, y <= min_y + epsilon ? min_y : y - epsilon};

            if (points_in_hull == 0) {
                first_x = x;
                rectangle[0] = p1;
                rectangle[1] = p2;
                upper.clear();
                lower.clear();
                upper.push_back(p1);
                lower.push_back(p2);
                upper_start = lower_start = 0;
                ++points_in_hull;
                return true;
            }

            if (points_in_hull == 1) {
                rectangle[2] = p2;
                rectangle[3] = p1;
                upper.push_back(p1);
                lower.push_back(p2);
                ++points_in_hull;
                return true;
            }

            auto slope1 = rectangle[2] - rectangle[0];
            auto slope2 = rectangle[3] - rectangle[1];
            bool outside_line1 = p1 - rectangle[2] < slope1;
            bool outside_line2 = p2 - rectangle[3] > slope2;

            if (outside_line1 || outside_line2) {
                points_in_hull = 0;
                return false;
            }

            if (p1 - rectangle[1] < slope2) {
                // Find extreme slope
                auto min = lower[lower_start] - p1;
                auto min_i = lower_start;
                for (auto i = lower_start + 1; i < lower.size(); i++) {
                    auto val = lower[i] - p1;
                    if (val > min)
                        break;
                    min = val;
                    min_i = i;
                }

                rectangle[1] = lower[min_i];
                rectangle[3] = p1;
                lower_start = min_i;

                // Hull update
                auto end = upper.size();
                for (; end >= upper_start + 2 && cross(upper[end - 2], upper[end - 1], p1) <= 0; --end)
                    continue;
                upper.resize(end);
                upper.push_back(p1);
            }

            if (p2 - rectangle[0] > slope1) {
                // Find extreme slope
                auto max = upper[upper_start] - p2;
                auto max_i = upper_start;
                for (auto i = upper_start + 1; i < upper.size(); i++) {
                    auto val = upper[i] - p2;
                    if (val < max)
                        break;
                    max = val;
                    max_i = i;
                }

                rectangle[0] = upper[max_i];
                rectangle[2] = p2;
                upper_start = max_i;

                // Hull update
                auto end = lower.size();
                for (; end >= lower_start + 2 && cross(lower[end - 2], lower[end - 1], p2) >= 0; --end)
                    continue;
                lower.resize(end);
                lower.push_back(p2);
            }

            ++points_in_hull;
            return true;
        }

        forceinline bool append_point(const X &x, const Y &y) {
//        if (points_in_hull > 0 && x <= last_x)
//            throw std::logic_error("Points must be increasing by x.");

            last_x = x;
            auto max_y = std::numeric_limits<Y>::max();
            auto min_y = std::numeric_limits<Y>::lowest();
            Point p1{x, y >= max_y - epsilon ? max_y : y + epsilon};
            Point p2{x, y <= min_y + epsilon ? min_y : y - epsilon};

            if (points_in_hull < 2) {
                buffqueue.push(p1);
                buffqueue.push(p2);
                if (buffqueue.size() == 4){
                    upper_start = lower_start = 0;
                    rectangle[0] =  buffqueue.front();
                    first_x = rectangle[0].x;
                    buffqueue.pop();
                    rectangle[1] = buffqueue.front();
                    buffqueue.pop();
                    rectangle[3] =  buffqueue.front();
                    buffqueue.pop();
                    rectangle[2] =  buffqueue.front();
                    buffqueue.pop();
                    upper.push_back(rectangle[0]);
                    lower.push_back(rectangle[1]);
                    upper.push_back(rectangle[3]);
                    lower.push_back(rectangle[2] );
                    points_in_hull ++;
                }
                else {
                    points_in_hull ++;
                }
//                else{
//                    cout<< "my Lord, i need You! please have pity on me!!"  << endl;
//                }

                return true;
            }

//            if (points_in_hull == 1) {
//                rectangle[2] = p2;
//                rectangle[3] = p1;
//                upper.push_back(p1);
//                lower.push_back(p2);
//                ++points_in_hull;
//                return true;
//            }

            auto slope1 = rectangle[2] - rectangle[0];
            auto slope2 = rectangle[3] - rectangle[1];
            bool outside_line1 = p1 - rectangle[2] < slope1;
            bool outside_line2 = p2 - rectangle[3] > slope2;

            if (outside_line1 || outside_line2) {
                points_in_hull = 0;
                vector<Point>().swap(lower);
                vector<Point>().swap(upper);
                lower_start = 0;
                upper_start = 0;
                return false;
            }

            if (p1 - rectangle[1] < slope2) {
                // Find extreme slope
                auto min = lower[lower_start] - p1;
                auto min_i = lower_start;
                for (auto i = lower_start + 1; i < lower.size(); i++) {
                    auto val = lower[i] - p1;
                    if (val > min)
                        break;
                    min = val;
                    min_i = i;
                }

                rectangle[1] = lower[min_i];
                rectangle[3] = p1;
                lower_start = min_i;

                // Hull update
                auto end = upper.size();
                for (; end >= upper_start + 2 && cross(upper[end - 2], upper[end - 1], p1) <= 0; --end)
                    continue;
                upper.resize(end);
                upper.push_back(p1);
            }

            if (p2 - rectangle[0] > slope1) {
                // Find extreme slope
                auto max = upper[upper_start] - p2;
                auto max_i = upper_start;
                for (auto i = upper_start + 1; i < upper.size(); i++) {
                    auto val = upper[i] - p2;
                    if (val < max)
                        break;
                    max = val;
                    max_i = i;
                }

                rectangle[0] = upper[max_i];
                rectangle[2] = p2;
                upper_start = max_i;

                // Hull update
                auto end = lower.size();
                for (; end >= lower_start + 2 && cross(lower[end - 2], lower[end - 1], p2) >= 0; --end)
                    continue;
                lower.resize(end);
                lower.push_back(p2);
            }

            ++points_in_hull;
            return true;
        }


        forceinline CanonicalSegment get_segment() {
            if (points_in_hull == 1)
                return CanonicalSegment(rectangle[0], rectangle[1], first_x);
            return CanonicalSegment(rectangle, first_x);
        }

        forceinline CanonicalSegment get_segment( X break_x) {
            if (points_in_hull == 1)
                return CanonicalSegment(rectangle[0], rectangle[1], first_x);
            return CanonicalSegment(rectangle, first_x,break_x);
        }

        forceinline CanonicalSegment get_segment( X break_x,std::vector<X> keys_vec) {
            if (points_in_hull == 1)
                return CanonicalSegment(rectangle[0], rectangle[1], first_x);
            return CanonicalSegment(rectangle, first_x,break_x,keys_vec);
        }


        void reset() {
            points_in_hull = 0;
            vector<Point>().swap(lower);
            vector<Point>().swap(upper);
            lower.clear();
            upper.clear();
        }
    };

    template<typename X, typename Y>
    class insertPWLF<X, Y>::CanonicalSegment {
        friend class insertPWLF;


    public:
        std::vector<X> slotkey;
        Point rectangle[4];
        X first;
        X last;

        CanonicalSegment( X first) : first(first) {};
        CanonicalSegment(const Point &p0, const Point &p1, X first) : rectangle{p0, p1, p0, p1}, first(first) {};

        CanonicalSegment(const Point &p0, const Point &p1, X first, X last) : rectangle{p0, p1, p0, p1}, first(first) ,last(last) {};

        CanonicalSegment(const Point (&rectangle)[4], X first)
                : rectangle{rectangle[0], rectangle[1], rectangle[2], rectangle[3]}, first(first) {};

        CanonicalSegment(const Point (&rectangle)[4], X first, X last)
                : rectangle{rectangle[0], rectangle[1], rectangle[2], rectangle[3]}, first(first), last(last) {};



        forceinline bool one_point() const {
            return rectangle[0].x == rectangle[2].x && rectangle[0].y == rectangle[2].y
                   && rectangle[1].x == rectangle[3].x && rectangle[1].y == rectangle[3].y;
        }

    public:

        CanonicalSegment() = default;

        X get_first_x() const { return first; }
        X get_last_x() const { return last; }

        forceinline std::pair<long double, long double> get_intersection() const {
            auto &p0 = rectangle[0];
            auto &p1 = rectangle[1];
            auto &p2 = rectangle[2];
            auto &p3 = rectangle[3];
            auto slope1 = p2 - p0;
            auto slope2 = p3 - p1;

            if (one_point() || slope1 == slope2)
                return {p0.x, p0.y};

            auto p0p1 = p1 - p0;
            auto a = slope1.dx * slope2.dy - slope1.dy * slope2.dx;
            auto b = (p0p1.dx * slope2.dy - p0p1.dy * slope2.dx) / static_cast<long double>(a);
            auto i_x = p0.x + b * slope1.dx;
            auto i_y = p0.y + b * slope1.dy;
            return {i_x, i_y};
        }

        forceinline std::pair<long double, SY> get_floating_point_segment(const X &origin) const {
            if (one_point())
                return {0, (rectangle[0].y + rectangle[1].y) / 2};

            if constexpr (std::is_integral_v<X> && std::is_integral_v<Y>) {
                auto slope = rectangle[3] - rectangle[1];
                auto intercept_n = slope.dy * (SX(origin) - rectangle[1].x);
                auto intercept_d = slope.dx;
                auto rounding_term = ((intercept_n < 0) ^ (intercept_d < 0) ? -1 : +1) * intercept_d / 2;
                auto intercept = (intercept_n + rounding_term) / intercept_d + rectangle[1].y;
                return {static_cast<long double>(slope), intercept};
            }

            auto[i_x, i_y] = get_intersection();
            auto[min_slope, max_slope] = get_slope_range();
            auto slope = (min_slope + max_slope) / 2.;
            auto intercept = i_y - (i_x - origin) * slope;
            return {slope, intercept};
        }

        forceinline std::pair<long double, long double> get_slope_range() const {
            if (one_point())
                return {0, 1};

            auto min_slope = static_cast<long double>(rectangle[2] - rectangle[0]);
            auto max_slope = static_cast<long double>(rectangle[3] - rectangle[1]);
            return {min_slope, max_slope};
        }
    };


    template<typename key_type,typename filmtype>
    std::pair<size_t,std::vector<key_type> > append_segmentation(size_t error,std::vector<key_type> keys,
                                                                 filmtype *filmada,unsigned int k){
        size_t c = 0;
        std::vector<key_type> startkeys;
         // innerpiece;
        for (size_t i = 0; i < keys.size(); ++i) {
            pair<key_type,int> p(keys[i],filmada->innerlevels[k]->nextpos++) ;   // i 为 pos
            ++(filmada->innerlevels[k]->pos);
            if (!filmada->innerlevels[k]->opt->add_point(p.first, p.second)) {  // 如果inner level  不满足error 了，那么再创建一个innerpiece

//                filmada->i_innerpiece = innerpiece;
                auto a = filmada->innerlevels[k]->opt->get_segment();
                filmada->innerlevels[k]->i_innerpiece->update(a);
                if (filmada->innerlevels[k]->pos > 2)
                    filmada->innerlevels[k]->innerpieces.pop_back();
                filmada->innerlevels[k]->innerpieces.emplace_back(filmada->innerlevels[k]->i_innerpiece);
                // 当前 innerpiece 不再满足，需要创建new inner piece 并判断该 innerlevel 的上一层level 是否需要更新
                typename filmtype:: Innerpiece* innerpiece = new typename filmtype:: Innerpiece();
                filmada->innerlevels[k]->i_innerpiece = innerpiece;
                // 首先在该层创建一个 new innerpiece， 更新该innerpiece，再递归向向上

                filmada->innerlevels[k]->pos = 0;
                filmada->innerlevels[k]->nextpos -= 2;
                delete filmada->innerlevels[k]->opt;
//                filmada->innerlevels[k]->opt.pop_back();
                insertPWLF<key_type, int> *inneropt = new insertPWLF<key_type, int>(error);
                filmada->innerlevels[k]->opt=inneropt;
                if (k==0){
//                    auto aaaaa = filmada->leaflevel.leafpieces.size()-2;
                    startkeys.emplace_back(filmada->leaflevel.leafpieces[filmada->leaflevel.leafpieces.size()-1]->startkey);
                }

                else{
//                    auto aaaaa = filmada->innerlevels[k-1]->innerpieces.size()-2;
                    startkeys.emplace_back(filmada->innerlevels[k-1]->innerpieces[filmada->innerlevels[k-1]->innerpieces.size()-2]->startkey);
                }
                startkeys.emplace_back(p.first);
                auto rr = append_segmentation(error,startkeys,filmada,k);

                ++c;
                if (filmada->innerlevels.back()->innerpieces.size() > 1)
                {
                    startkeys.clear();
                    startkeys.emplace_back(a.first);
                    startkeys.emplace_back(p.first);
                    insertPWLF<key_type, int> *inneropt = new insertPWLF<key_type, int>(error);
//
                    typename filmtype::Innerlevel *innerlevel = new typename filmtype::Innerlevel;
                    filmada->innerlevels.emplace_back(innerlevel);
                    filmada->innerlevels.back()->opt = inneropt ;
                    auto rr = append_segmentation(error,startkeys,filmada,k+1);
//                    cout<< "Jesus, i need You !"  << endl;
                    startkeys.emplace_back(a.first);
                    return std::pair<size_t,std::vector<key_type> >(++c,startkeys) ;
                }
                else if (filmada->innerlevels.size() > 1 && (k != filmada->innerlevels.size()-1) )   // 如上为 由于创建了new innner piece，导致了new innerlevel，如下为，虽然创建了new innerpiece，但只需要更新上层的innner level
                {
                    // 更新上层的最后一个inner piece
                    startkeys.pop_back();
                    auto rr = append_segmentation(error,startkeys,filmada,k+1);
                    startkeys.clear();
//                    cout << "thank You, my Lord! i need You!"<<endl;
                    startkeys.emplace_back(a.first);
                    return std::pair<size_t,std::vector<key_type> >(++c,startkeys) ;
                }
                else if (filmada->innerlevels.size() > 1 && (k == filmada->innerlevels.size()-1) )   // 如上为 由于创建了new innner piece，导致了new innerlevel，如下为，虽然创建了new innerpiece，但只需要更新上层的innner level
                {
                    cout << "thank You, my Lord! i need You!"<<endl;
                    startkeys.clear();
                    return std::pair<size_t,std::vector<key_type> >(++c,startkeys) ;
                }
                else{
                    startkeys.clear();
                }
                a = filmada->innerlevels[k]->opt->get_segment();
                filmada->innerlevels[k]->i_innerpiece->update(a);
                if (filmada->innerlevels[k]->pos > 2)
                    filmada->innerlevels[k]->innerpieces.pop_back();
                filmada->innerlevels[k]->innerpieces.emplace_back(filmada->innerlevels[k]->i_innerpiece);

                startkeys.emplace_back(a.first);
                return std::pair<size_t,std::vector<key_type> >(++c,startkeys) ;
            }

        }

        auto a = filmada->innerlevels[k]->opt->get_segment();
        filmada->innerlevels[k]->i_innerpiece->update(a);
        if (filmada->innerlevels[k]->pos > 2)
            filmada->innerlevels[k]->innerpieces.pop_back();
        filmada->innerlevels[k]->innerpieces.emplace_back(filmada->innerlevels[k]->i_innerpiece);

        startkeys.emplace_back(a.first);
        return std::pair<size_t,std::vector<key_type> >(++c,startkeys) ;

    }

//
    template<typename key_type, typename leaf_type,typename filmtype,typename filmadalrutype>
    std::pair<size_t,std::vector<key_type> > append_segmentation(size_t error,std::vector<key_type> keys,key_type* payload,
                                                                 filmtype *filmada, unsigned long int &inkeynum,leaf_type* m_tailleaf,filmadalrutype *interchain){


        size_t c = 0;
        size_t start = 0;
        std::vector<key_type> startkeys;
        unsigned int pos = 0;
        pair<key_type,unsigned int> p(keys[0],pos++) ;
//        if (keys[0] == 62129282194 ){
//            cout<< "Jesus, please come!, 62129282194" <<endl;
//        }
//        int valuesize = filmada->valuesize;
        inkeynum +=1;


        if (m_tailleaf == NULL)  // 初始化
        {
            leaf_type *cur_leaf = new leaf_type;
            cur_leaf->slotkey.reserve(8192*4);
            cur_leaf->slotdata.reserve(8192*4);
            cur_leaf->locbitmap.reserve(8192*4);
            cur_leaf->delbitmap.reserve(8192*4);
            filmada->m_tailleaf = cur_leaf;
            filmada->leaflevel.leafpieces.emplace_back( filmada->m_tailleaf);
//            m_tailleaf = cur_leaf;
            insertPWLF<key_type, int> *leafopt = new insertPWLF<key_type, int> (error);
            filmada->leaflevel.opt = leafopt ;

        }
//        auto xx = filmada->leaflevel.second[0]->buffqueue.size();
        filmada->leaflevel.opt->append_point(p.first, p.second);
//        key_type* value = new key_type[valuesize];
        filmada->m_tailleaf->slotdata.emplace_back( filmada->m_tailleaf->intrachain.put(p.second,payload));
        filmada->m_tailleaf->slotkey.emplace_back(p.first);
        filmada->m_tailleaf->locbitmap.emplace_back(true);
        filmada->m_tailleaf->delbitmap.emplace_back(false);
        for (inkeynum; inkeynum < keys.size(); ++inkeynum) {
//            if (keys[inkeynum] == 142819719223789952 ){
//                cout<< "Jesus, please come!, 142819719223789952" <<endl;
//            }
            pair<key_type,unsigned int> next_p(keys[inkeynum],pos++) ;

            if (inkeynum != start && next_p.first == p.first)
                continue;
            p = next_p;
            if (filmada ->leaflevel.opt->append_point(p.first, p.second)) {
                filmada->m_tailleaf->slotkey.emplace_back(p.first);
//                key_type* value = new key_type[valuesize];
                filmada->m_tailleaf->slotdata.emplace_back( filmada->m_tailleaf->intrachain.put(p.second,payload));
                filmada->m_tailleaf->locbitmap.emplace_back( true);
                filmada->m_tailleaf->delbitmap.emplace_back( false);
            }
            else
            {
                start = inkeynum;
                auto a = filmada ->leaflevel.opt->get_segment(keys[--inkeynum]); // 将生成的new leaf piece插入到leaflevel 中
                filmada->m_tailleaf->update(a);
//                filmada ->leaflevel.leafpieces.emplace_back( filmada->m_tailleaf);
                interchain->put(filmada->m_tailleaf->startkey,filmada->m_tailleaf);
//                if (filmada->m_tailleaf->startkey == 12367821099355)
//                    COUT_THIS('thank You, my lovely Lord');
                auto buffersize = ceil(filmada->m_tailleaf->slotkey.size()*filmada->buffer_ratio)+1;   //0.429, 0.25,0.667
                filmada->m_tailleaf->buffer_piece = new typename filmtype:: sortPiece(buffersize);

                if (buffersize < filmada->minbsize)
                    filmada->minbsize = buffersize;
                else if (buffersize > filmada->maxbsize)
                    filmada->maxbsize = buffersize;
                filmada->m_tailleaf->buffer_piece->slope = buffersize;
                filmada->m_tailleaf->buffer_piece->slotkey.reserve(buffersize);
                filmada->m_tailleaf->buffer_piece->slotdata.reserve(buffersize);
                filmada->m_tailleaf->buffer_piece->locbitmap.reserve(buffersize);
                // 这里是初始化 parent piece 的first key 和 last key
                if (filmada->innerlevels.size() == 0){
                    startkeys.emplace_back( filmada->m_tailleaf->startkey);
                    startkeys.emplace_back( p.first);
//                    typename filmtype:: Innerpiece* innerpiece = new typename filmtype:: Innerpiece() ;// 创建parent piece
//                    filmada->i_innerpiece = innerpiece;
                    insertPWLF<key_type, int> *inneropt = new insertPWLF<key_type, int>(error);
                    typename filmtype::Innerlevel *innerlevel = new typename filmtype::Innerlevel;
                    filmada->innerlevels.emplace_back(innerlevel);
                    filmada->innerlevels[0]->opt = inneropt ;
//                    cout << " my Lord, Jesus, please have pity on me"<< endl;
                }
                else{
                    // 从innerlevel 的最底层到root 层，判断是否需要更新
                    startkeys.emplace_back( p.first);
//                    cout << "my lovely Lord, i trust in You!" << endl;
                }
                auto rr = append_segmentation(error,startkeys,filmada,0);
//                startkeys.clear();
                vector<key_type>().swap(startkeys);
//                cout << "Jesus, i need You!!"<< endl;
                pos = 0;

                filmada->m_tailleaf = new leaf_type;
                filmada->m_tailleaf->slotkey.reserve(8192*4);
                filmada->m_tailleaf->slotdata.reserve(8192*4);
                filmada->m_tailleaf->locbitmap.reserve(8192*4);
                filmada->m_tailleaf->delbitmap.reserve(8192*4);
                filmada ->leaflevel.leafpieces.emplace_back( filmada->m_tailleaf);
                ++c;
            }
        }
        auto a = filmada ->leaflevel.opt->get_segment( keys.back());
        if (filmada ->leaflevel.opt->points_in_hull < 2){
            a.first = p.first;
            a.last = p.first;
        }
        filmada->m_tailleaf->update(a);
        interchain->put(filmada->m_tailleaf->startkey,filmada->m_tailleaf);
        auto buffersize = (filmada->maxbsize+filmada->maxbsize)/2;
        filmada->m_tailleaf->buffer_piece = new typename filmtype:: sortPiece(buffersize);
        filmada->m_tailleaf->buffer_piece->slope = buffersize;
        filmada->m_tailleaf->buffer_piece->slotkey.reserve(buffersize);
        filmada->m_tailleaf->buffer_piece->slotdata.reserve(buffersize);
        filmada->m_tailleaf->buffer_piece->locbitmap.reserve(buffersize);
        return std::pair<size_t,vector<key_type>> (++c,startkeys);
    }

    // single piece retrain，use the idea of link. when the model of new models >= 2, inserted to the buffer piece
    template<typename key_type, typename filmtype,typename filmadalrutype>
    typename filmtype::sortPiece* make_segmentation(filmtype *filmindex,size_t error,std::vector<key_type> retrainkeys,
                                                    std::vector<bool> retrainflags, std::vector<void*> retraindata,
                                                    filmadalrutype *interchain,int leafslot){

//        size_t c = 0;
        size_t start = 0;
        unsigned int pos = 0;
        vector<typename filmtype::Leafpiece*> newleafs;
        typename filmtype::Leafpiece* curleaf = new typename filmtype::Leafpiece;
        insertPWLF<key_type, int> *leafopt = new insertPWLF<key_type, int> (error);
//        if (keys[0] == 17386164347307 ){
//            cout<< "i need You, my Lord, rebuild the model~~~" <<endl;
//        }

        for (int i = 0; i < retrainkeys.size(); i++) {
//            if (keys[i] == 17386164347307 ) { //  15945001399101
//                cout << "my Lord, i need You! thank You, for ever and ever ~~" << endl;
//            }
            pair<key_type,unsigned int> p(retrainkeys[i],pos++);
            if (leafopt->append_point(p.first, p.second)) {
                curleaf->slotkey.emplace_back(retrainkeys[i]);
                if (retrainflags[i] == true){
                    curleaf->locbitmap.emplace_back(true);
                    auto finddata = (adalru::Node<lruOff_type, key_type* >*) retraindata[i];
                    curleaf->slotdata.emplace_back( curleaf->intrachain.put(p.second,finddata->value));
                }
                else{
                    curleaf->slotdata.emplace_back(retraindata[i]);
                    curleaf->locbitmap.emplace_back(false);
                }
            }
            else
            {
                start = i;
                auto a = leafopt->get_segment(retrainkeys[--i]); // 将生成的new leaf piece插入到leaflevel 中
                curleaf->update(a);
                newleafs.emplace_back(curleaf);
                if (curleaf->intrachain.size >0)
                    interchain->put(curleaf->startkey,curleaf);

                if (start >=  retrainkeys.size() - 4){
//                    cout << "i need you, my Lord, the next leaf is too small, so the remaining key should be inserted into the buffer piece"<< endl;
                    auto buffersize = ceil(newleafs[0]->slotkey.size()*filmindex->buffer_ratio)+4;
                    typename filmtype::sortPiece* bufferpiece = new typename filmtype:: sortPiece(buffersize);
                    curleaf->buffer_piece = bufferpiece;
                    i++;
                    for (; i < retrainkeys.size(); i++){
                        auto key = retrainkeys[i];
//                        if (key == 15143644349635  )
//                            cout << "i need You, my lovely Lord, rebuild 17386164347307" << endl;
                        auto predicted_pos = (*newleafs[0])(key)*filmindex->buffer_ratio;
                        if (retrainflags[i] == true){
                            auto finddata = (adalru::Node<lruOff_type, key_type* >*) retraindata[i];
                            bufferpiece->model_insert(key,true,finddata->value,predicted_pos);
                        }
                        else
                            bufferpiece->model_insert(key,false,retraindata[i],predicted_pos);
                    }
//                    bufferpiece->data_capacity_++;
//                    if (bufferpiece->num_keys_ > 0){
//                        interchain->put(bufferpiece->ALEX_DATA_NODE_KEY_AT(0),bufferpiece);
//                    }
                }
                else{
                    auto buffersize = ceil(newleafs[0]->slotkey.size()*filmindex->buffer_ratio)+1;
                    typename filmtype::sortPiece* bufferpiece = new typename filmtype:: sortPiece(buffersize);
                    curleaf->buffer_piece = bufferpiece;
                    pos = 0;
                    delete leafopt;
                    leafopt = NULL;
                    leafopt = new insertPWLF<key_type, int> (error);
                    curleaf = new typename filmtype::Leafpiece;
//                    ++c;
                }
            }
        }

        if (curleaf->buffer_piece == NULL){
            auto a = leafopt->get_segment(retrainkeys.back()); // 将生成的new leaf piece插入到leaflevel 中
            curleaf->update(a);
            delete leafopt;
            leafopt = NULL;
            if (curleaf->intrachain.size >0)
                interchain->put(curleaf->slotkey[0],curleaf);
//            else
//                cout << "Jesus, i need You, curleaf->intrachain.size == 0" << endl;
            newleafs.emplace_back(curleaf);
            auto buffersize = ceil(newleafs[0]->slotkey.size()*filmindex->buffer_ratio)+1;
            curleaf->buffer_piece = new typename filmtype:: sortPiece(buffersize);
//            auto buffersize = ceil(curleaf->slotkey.size()*0.667);

//            curleaf->buffer_piece->slotkey.reserve(buffersize);
//            curleaf->buffer_piece->slotdata.reserve(buffersize);
//            curleaf->buffer_piece->locbitmap.reserve(buffersize);

//            ++c;
        }

        // new leafs 要替换原来find leaf 的位置
        // 1.将leafslot 位置的leaf 替换为newleaf[0]
        filmindex->leaflevel.leafpieces[leafslot] = newleafs[0];
//        filmtype::delete_leaf(leaf); // 释放original leaf
//        auto original1 = filmindex->leaflevel.leafpieces[leafslot];
//        if (newleafs[0]->startkey == 15143643742681 ) { //  15945001399101
//            cout << "my Lord, i need You! thank You, for ever and ever ~~" << endl;
//        }
        // 2. if newleafs has more than one leaf, link to the buffer piece
        auto c = newleafs.size();
        if (c>1){
            typename filmtype::Leafpiece* curleaf = newleafs[0];
            for (int li = 1;li<c;li++){
//                if (curleaf->startkey == 15143643742681 ) { //  15945001399101
//                    cout << "my Lord, i need You! thank You, for ever and ever ~~" << endl;
//                }
                curleaf->buffer_piece = newleafs[li];
                curleaf = newleafs[li];
            }
        }
//        if (curleaf->buffer_piece->buffer_piece!=NULL)
//            cout << "Jeses, i need You! c = 0"<< endl;


        //
        return  (typename filmtype:: sortPiece*)curleaf->buffer_piece;

    }

    // during inner rebuild
    // single piece retrain，use the idea of link. when the model of new models >= 2, inserted to the buffer piece
    template<typename key_type, typename filmtype,typename filmadalrutype>
    typename filmtype::sortPiece* make_segmentation(filmtype *filmindex,size_t error,std::vector<key_type> retrainkeys,
                                                    std::vector<bool> retrainflags, std::vector<void*> retraindata,
                                                    filmadalrutype *interchain,vector<typename filmtype::Leafpiece*>& linkpieces, int li){

//        size_t c = 0;
        size_t start = 0;
        unsigned int pos = 0;
        vector<typename filmtype::Leafpiece*> newleafs;
        typename filmtype::Leafpiece* curleaf = new typename filmtype::Leafpiece;
        insertPWLF<key_type, int> *leafopt = new insertPWLF<key_type, int> (error);
//        if (keys[0] == 17386164347307 ){
//            cout<< "i need You, my Lord, rebuild the model~~~" <<endl;
//        }

        for (int i = 0; i < retrainkeys.size(); i++) {
//            if (keys[i] == 17386164347307 ) { //  15945001399101
//                cout << "my Lord, i need You! thank You, for ever and ever ~~" << endl;
//            }
            pair<key_type,unsigned int> p(retrainkeys[i],pos++);
            if (leafopt->append_point(p.first, p.second)) {
                curleaf->slotkey.emplace_back(retrainkeys[i]);
                if (retrainflags[i] == true){
                    curleaf->locbitmap.emplace_back(true);
                    auto finddata = (adalru::Node<lruOff_type, key_type* >*) retraindata[i];
                    curleaf->slotdata.emplace_back( curleaf->intrachain.put(p.second,finddata->value));
                }
                else{
                    curleaf->slotdata.emplace_back(retraindata[i]);
                    curleaf->locbitmap.emplace_back(false);
                }
            }
            else
            {
                start = i;
                auto a = leafopt->get_segment(retrainkeys[--i]); // 将生成的new leaf piece插入到leaflevel 中
                curleaf->update(a);
                newleafs.emplace_back(curleaf);
                if (curleaf->intrachain.size >0)
                    interchain->put(curleaf->startkey,curleaf);

                if (start >=  retrainkeys.size() - 4){
//                    cout << "i need you, my Lord, the next leaf is too small, so the remaining key should be inserted into the buffer piece"<< endl;
                    auto buffersize = ceil(newleafs[0]->slotkey.size()*filmindex->buffer_ratio)+4;
                    typename filmtype::sortPiece* bufferpiece = new typename filmtype:: sortPiece(buffersize);
                    curleaf->buffer_piece = bufferpiece;
                    i++;
                    for (; i < retrainkeys.size(); i++){
                        auto key = retrainkeys[i];
//                        if (key == 15143644349635  )
//                            cout << "i need You, my lovely Lord, rebuild 17386164347307" << endl;
                        auto predicted_pos = (*newleafs[0])(key)*filmindex->buffer_ratio;
                        if (retrainflags[i] == true){
                            auto finddata = (adalru::Node<lruOff_type, key_type* >*) retraindata[i];
                            bufferpiece->model_insert(key,true,finddata->value,predicted_pos);
                        }
                        else
                            bufferpiece->model_insert(key,false,retraindata[i],predicted_pos);
                    }
//                    bufferpiece->data_capacity_++;
//                    if (bufferpiece->num_keys_ > 0){
//                        interchain->put(bufferpiece->ALEX_DATA_NODE_KEY_AT(0),bufferpiece);
//                    }
                }
                else{
                    auto buffersize = ceil(newleafs[0]->slotkey.size()*filmindex->buffer_ratio)+1;
                    typename filmtype::sortPiece* bufferpiece = new typename filmtype:: sortPiece(buffersize);
                    curleaf->buffer_piece = bufferpiece;
                    pos = 0;
                    delete leafopt;
                    leafopt = NULL;
                    leafopt = new insertPWLF<key_type, int> (error);
                    curleaf = new typename filmtype::Leafpiece;
//                    ++c;
                }
            }
        }

        if (curleaf->buffer_piece == NULL){
            auto a = leafopt->get_segment(retrainkeys.back()); // 将生成的new leaf piece插入到leaflevel 中
            curleaf->update(a);
            delete leafopt;
            leafopt = NULL;
            if (curleaf->intrachain.size >0)
                interchain->put(curleaf->slotkey[0],curleaf);
//            else
//                cout << "Jesus, i need You, curleaf->intrachain.size == 0" << endl;
            newleafs.emplace_back(curleaf);
            auto buffersize = ceil(newleafs[0]->slotkey.size()*filmindex->buffer_ratio)+1;
            curleaf->buffer_piece = new typename filmtype:: sortPiece(buffersize);
//            auto buffersize = ceil(curleaf->slotkey.size()*0.667);

//            curleaf->buffer_piece->slotkey.reserve(buffersize);
//            curleaf->buffer_piece->slotdata.reserve(buffersize);
//            curleaf->buffer_piece->locbitmap.reserve(buffersize);

//            ++c;
        }

        // new leafs 要替换原来find leaf 的位置
        // 1.将leafslot 位置的leaf 替换为newleaf[0]

        linkpieces[li] = newleafs[0];
//        filmtype::delete_leaf(leaf); // 释放original leaf
//        auto original1 = filmindex->leaflevel.leafpieces[leafslot];
//        if (newleafs[0]->startkey == 15143643742681 ) { //  15945001399101
//            cout << "my Lord, i need You! thank You, for ever and ever ~~" << endl;
//        }
        // 2. if newleafs has more than one leaf, link to the buffer piece
        auto c = newleafs.size();
        if (c>1){
            typename filmtype::Leafpiece* curleaf = newleafs[0];
            for (int li = 1;li<c;li++){
//                if (curleaf->startkey == 15143643742681 ) { //  15945001399101
//                    cout << "my Lord, i need You! thank You, for ever and ever ~~" << endl;
//                }
                curleaf->buffer_piece = newleafs[li];
                curleaf = newleafs[li];
            }
        }
//        if (curleaf->buffer_piece->buffer_piece!=NULL)
//            cout << "Jeses, i need You! c = 0"<< endl;


        //
        return  (typename filmtype:: sortPiece*)curleaf->buffer_piece;

    }

    // build the parent of leaf level, unlink the linked leaf pieces
    template< typename filmtype ,typename key_type,typename leaf_level_type,typename innerlevel_type, typename lru_type>
    std::pair<unsigned int ,vector<key_type>>  make_segmentation(filmtype* filmindex,unsigned short error,leaf_level_type *leaflevel, innerlevel_type *innerlevel, key_type firstkey, lru_type interchain) {

        unsigned int c = 0;
        unsigned int start = 0;
        vector<key_type> startkeys;
        pair<key_type,int> p;
        insertPWLF<key_type, int> *inneropt = new insertPWLF<key_type, int>(error);
        innerlevel->opt = inneropt ;
        vector<typename filmtype::Leafpiece*> leafpieces;

        for (unsigned int i = 0; i < leaflevel->leafpieces.size(); ++i) {
            ++(innerlevel->pos);
            // 依次添加 i 位置上的所有learned pieces （linked learned pieces）
            vector<typename filmtype::Leafpiece*> linkpieces;

            auto curleaf = leaflevel->leafpieces[i];
            if (curleaf->startkey == 16182388766135)
                cout << "Jesus, i need You, curleaf->startkey == 16182388766135" << endl;
            linkpieces.emplace_back(curleaf);
            if (curleaf->buffer_piece->buffer_piece) // 如果有link piece
            {
                if (curleaf->startkey == 16182388766135)
                    cout << "Jesus, i need You, curleaf->startkey == 16182388766135" << endl;
                while (curleaf->buffer_piece->buffer_piece){
                    curleaf = curleaf->buffer_piece;
                    linkpieces.emplace_back(curleaf);
                }
                // 找到该link 上的 最终的buffer piece
                typename filmtype::sortPiece* bufferpiece = (typename filmtype::sortPiece*)curleaf->buffer_piece;
                // 将buffer piece 中的 每个key 从新插入到linkpieces 中
                int li = 0;
//                if (bufferpiece->num_keys_>0 && bufferpiece->ALEX_DATA_NODE_KEY_AT(0)>linkpieces.back()->startkey){
//                    // 只需要将除 最后一个 leaf 以外的所有 buffer 都 new 一个 sort piece
//                    for (int pi = 0; pi < linkpieces.size()-1; pi++){
//                        auto buffersize = ceil(linkpieces[pi]->slotkey.size()*filmindex->buffer_ratio);
//                        linkpieces[pi]->buffer_piece = new typename filmtype::sortPiece(buffersize);
////                        if (buffersize < filmindex->minbsize)
////                            filmindex->minbsize = buffersize;
////                        else if (buffersize > filmindex->maxbsize)
////                            filmindex->maxbsize = buffersize;
//                    }
//                }

                // 每一个 leaf 的 buffer 都 new 一个 sort piece
                for (auto piece: linkpieces){
                    auto buffersize = ceil(piece->slotkey.size()*filmindex->buffer_ratio)+4;
                    piece->buffer_piece = new typename filmtype::sortPiece(buffersize);
//                    piece->buffer_piece->slope = buffersize;
//                        if (buffersize < filmindex->minbsize)
//                            filmindex->minbsize = buffersize;
//                        else if (buffersize > filmindex->maxbsize)
//                            filmindex->maxbsize = buffersize;
                }


                if (bufferpiece->num_keys_>0)
                    interchain->remove(bufferpiece->ALEX_DATA_NODE_KEY_AT(0));
                for (int qi = 0; qi < bufferpiece->data_capacity_; qi ++){
                    if (bufferpiece->check_exists(qi)){
                        // 确定要插入哪个piece
                        auto bufferkey = bufferpiece->ALEX_DATA_NODE_KEY_AT(qi);
//                        if (bufferkey == 16182566776953)
//                            COUT_THIS("thank You, my Lord");
                        for (; li<linkpieces.size();li++){
                            if (bufferkey <= linkpieces[li]->startkey)
                                break;}
                        li--;
                        typename filmtype::sortPiece* buffer = (typename filmtype::sortPiece*) (linkpieces[li]->buffer_piece);
//                        buffer->data_capacity_ +=1;
                        auto predicted_pos = (*linkpieces[li])(bufferkey)*filmindex->buffer_ratio;
                        buffer->model_insert(bufferkey,bufferpiece->ALEX_DATA_NODE_FLAG_AT(qi),bufferpiece->ALEX_DATA_NODE_PAYLOAD_AT(qi),predicted_pos);
                        if (buffer->num_keys_ == buffer->data_capacity_){
                            std::vector<key_type> retrain_keys;
                            std::vector<void*> retrain_datas;
                            std::vector<bool> retrain_flags;
                            linkpieces[li]->mlearnedMergeSort(error, interchain,retrain_keys,retrain_datas,retrain_flags);

                            // retrain model, get the new slope and intercept
                            bufferpiece = internal::make_segmentation( filmindex,error,retrain_keys,retrain_flags,retrain_datas,interchain,linkpieces,li) ;
//                this->newmodels += (size_buffer.first- link_num);
                            // 执行插入

                            predicted_pos = (*linkpieces[li])(bufferkey) * filmindex->buffer_ratio;  // normalize to buffer size;
                            bufferpiece->model_insert(bufferkey,bufferpiece->ALEX_DATA_NODE_FLAG_AT(qi),bufferpiece->ALEX_DATA_NODE_PAYLOAD_AT(qi),predicted_pos);

                        }

                    }
                }
                // 释放/删除buffer piece
//                delete bufferpiece;

                for (auto leaf:linkpieces){
                    leafpieces.emplace_back(leaf);
                    typename filmtype::sortPiece* buffer = (typename filmtype::sortPiece*) (leaf->buffer_piece);
                    // 将 buffer piece 插入到 global chain 中
                    if (buffer->num_keys_ > 0)
                        interchain->put(buffer->ALEX_DATA_NODE_KEY_AT(0), buffer);
                    p = pair<key_type,int> (leaf->startkey,innerlevel->nextpos++) ;
                    if (!innerlevel->opt->add_point(p.first, p.second)) {
//                        --i;
                        auto a = innerlevel->opt->get_segment(); // 将生成的new leaf piece插入到leaflevel 中
                        innerlevel->i_innerpiece->update(a);
                        innerlevel->innerpieces.emplace_back(innerlevel->i_innerpiece);
                        // 当前 innerpiece 不再满足，需要创建new inner piece 并判断该 innerlevel 的上一层level 是否需要更新
                        typename filmtype:: Innerpiece* innerpiece = new typename filmtype:: Innerpiece();
                        innerlevel->i_innerpiece = innerpiece;
                        // 首先在该层创建一个 new innerpiece， 更新该innerpiece，再递归向向上

                        innerlevel->pos = 0;
                        innerlevel->nextpos -= 2;
                        delete innerlevel->opt;
//                filmada->innerlevels[k]->opt.pop_back();
                        insertPWLF<key_type, int> *inneropt = new insertPWLF<key_type, int>(error);
                        innerlevel->opt=inneropt;
                        startkeys.emplace_back(a.first);
                        ++c;
                    }
                }
            }
            else{
                leafpieces.emplace_back(leaflevel->leafpieces[i]);
                p = pair<key_type,int> (leaflevel->leafpieces[i]->startkey,innerlevel->nextpos++) ;
                if (!innerlevel->opt->add_point(p.first, p.second)) {

//                    --i;
                    auto a = innerlevel->opt->get_segment(); // 将生成的new leaf piece插入到leaflevel 中
                    innerlevel->i_innerpiece->update(a);
                    innerlevel->innerpieces.emplace_back(innerlevel->i_innerpiece);
                    // 当前 innerpiece 不再满足，需要创建new inner piece 并判断该 innerlevel 的上一层level 是否需要更新
                    typename filmtype:: Innerpiece* innerpiece = new typename filmtype:: Innerpiece();
                    innerlevel->i_innerpiece = innerpiece;
                    // 首先在该层创建一个 new innerpiece， 更新该innerpiece，再递归向向上

                    innerlevel->pos = 0;
                    innerlevel->nextpos -= 2;
                    delete innerlevel->opt;
//                filmada->innerlevels[k]->opt.pop_back();
                    insertPWLF<key_type, int> *inneropt = new insertPWLF<key_type, int>(error);
                    innerlevel->opt=inneropt;
                    startkeys.emplace_back(a.first);
                    ++c;
                }
            }

        }

        auto a = innerlevel->opt->get_segment();
        innerlevel->i_innerpiece->update(a);
        innerlevel->innerpieces.emplace_back(innerlevel->i_innerpiece);
        vector<typename filmtype::Leafpiece*>().swap(leaflevel->leafpieces);
        leaflevel->leafpieces = leafpieces;
        startkeys.emplace_back(a.first);

        return std::pair<unsigned int,vector<key_type> >(++c,startkeys) ;
    }



    // build the parent of leaf level
    template< typename filmtype ,typename key_type,typename leaf_level_type,typename innerlevel_type>
    std::pair<unsigned int ,vector<key_type>>  make_segmentation(filmtype* filmindex,unsigned short error,leaf_level_type *leaflevel, innerlevel_type *innerlevel, key_type firstkey) {

        unsigned int c = 0;
        unsigned int start = 0;
        vector<key_type> startkeys;
        pair<key_type,int> p(firstkey,innerlevel->nextpos++) ;
        insertPWLF<key_type, int> *inneropt = new insertPWLF<key_type, int>(error);
        innerlevel->opt = inneropt ;
        innerlevel->opt->add_point(p.first, p.second);

        for (unsigned int i = 1; i < leaflevel->leafpieces.size(); ++i) {
            ++(innerlevel->pos);
            pair<key_type,int> next_p(leaflevel->leafpieces[i]->startkey,innerlevel->nextpos++) ;
//            if (i != start && next_p.first == p.first)
//                continue;
            p = next_p;
            if (!innerlevel->opt->add_point(p.first, p.second)) {

                start = i;
                --i;
                auto a = innerlevel->opt->get_segment(); // 将生成的new leaf piece插入到leaflevel 中
                innerlevel->i_innerpiece->update(a);
                innerlevel->innerpieces.emplace_back(innerlevel->i_innerpiece);
                // 当前 innerpiece 不再满足，需要创建new inner piece 并判断该 innerlevel 的上一层level 是否需要更新
                typename filmtype:: Innerpiece* innerpiece = new typename filmtype:: Innerpiece();
                innerlevel->i_innerpiece = innerpiece;
                // 首先在该层创建一个 new innerpiece， 更新该innerpiece，再递归向向上

                innerlevel->pos = 0;
                innerlevel->nextpos -= 2;
                delete innerlevel->opt;
//                filmada->innerlevels[k]->opt.pop_back();
                insertPWLF<key_type, int> *inneropt = new insertPWLF<key_type, int>(error);
                innerlevel->opt=inneropt;
                startkeys.emplace_back(a.first);

                ++c;
            }
        }

        auto a = innerlevel->opt->get_segment();
        innerlevel->i_innerpiece->update(a);
        innerlevel->innerpieces.emplace_back(innerlevel->i_innerpiece);
        startkeys.emplace_back(a.first);

        return std::pair<unsigned int,vector<key_type> >(++c,startkeys) ;
    }


    // internal level retrain
    template<typename filmtype, typename key_type,typename innerlevel_type>
    std::pair<unsigned int,vector<key_type> > make_segmentation(filmtype *filmindex,size_t error,std::vector<key_type> keys,innerlevel_type *innerlevel){

        unsigned int c = 0;
        unsigned int start = 0;
        vector<key_type> startkeys;
        pair<key_type,int> p(keys[0],innerlevel->nextpos++) ;
        insertPWLF<key_type, int> *inneropt = new insertPWLF<key_type, int>(error);
        innerlevel->opt = inneropt ;
        innerlevel->opt->add_point(p.first, p.second);

        for (unsigned int i = 1; i < keys.size(); ++i) {
            ++(innerlevel->pos);
            pair<key_type,int> next_p(keys[i],innerlevel->nextpos++) ;
            if (i != start && next_p.first == p.first)
                continue;
            p = next_p;
            if (!innerlevel->opt->add_point(p.first, p.second)) {

                start = i;
                --i;
                auto a = innerlevel->opt->get_segment(); // 将生成的new leaf piece插入到leaflevel 中
                innerlevel->i_innerpiece->update(a);
                innerlevel->innerpieces.emplace_back(innerlevel->i_innerpiece);
                startkeys.emplace_back(a.first);
                // 当前 innerpiece 不再满足，需要创建new inner piece 并判断该 innerlevel 的上一层level 是否需要更新
                typename filmtype:: Innerpiece* innerpiece = new typename filmtype:: Innerpiece();
                innerlevel->i_innerpiece = innerpiece;
                // 首先在该层创建一个 new innerpiece， 更新该innerpiece，再递归向向上

                innerlevel->pos = 0;
                innerlevel->nextpos -= 2;
                delete innerlevel->opt;
//                filmada->innerlevels[k]->opt.pop_back();
                insertPWLF<key_type, int> *inneropt = new insertPWLF<key_type, int>(error);
                innerlevel->opt=inneropt;
                ++c;
            }
        }

        auto a = innerlevel->opt->get_segment();
        innerlevel->i_innerpiece->update(a);
        innerlevel->innerpieces.emplace_back(innerlevel->i_innerpiece);
        startkeys.emplace_back(a.first);

        vector<key_type>().swap(keys);
        return std::pair<unsigned int,vector<key_type> >(++c,startkeys) ;
    }



}




#endif //FILMINSERT_PWLF_H
