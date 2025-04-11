//
// Created by chaochao on 2021/12/23.
//

#ifndef EXPERIMENTCC12_FILMADALRU_H
#define EXPERIMENTCC12_FILMADALRU_H
#define forceinline inline __attribute__((__always_inline__))
#include <unordered_map>

#define MAX_INT    (((unsigned int)(-1))>>1)
// define the node in doubly-linked list
using namespace std;
namespace adalru {
    template<class Key, class Value>
    struct Node {
        Key key;
        Value value;
        Node *prev;
        Node *next;

        Node(Key k, Value v) : key(k), value(v), prev(nullptr), next(nullptr) {};

        Node() : key(), value(), prev(nullptr), next(nullptr) {};
        ~Node(){
            prev = nullptr;
            next = nullptr;
        }
//    Node(): prev(nullptr), next(nullptr){ };
    };

    template<class Key>
    struct bufNode {
        Key key;
        bufNode *prev;
        bufNode *next;

        bufNode(Key k) : key(k), prev(nullptr), next(nullptr) {};

        bufNode() : key(), prev(nullptr), next(nullptr) {};
        ~bufNode(){
            prev = nullptr;
            next = nullptr;
        }
//    Node(): prev(nullptr), next(nullptr){ };
    };

    template<class Key, class Value, class mapvalue>   //mapvalue 是指向在 lru 中的node 的指针  在hash map 中使用
    class hashLRU {

    public:
        int size = 0;
        int capacity;
        std::unordered_map<Key, mapvalue> map;
        Node<Key, Value> *head;
        Node<Key, Value> *tail;

        hashLRU(int def_capcity) {
            capacity = def_capcity;
            head = new Node<Key, Value>;
            tail = new Node<Key, Value>;
            head->next = tail;
            tail->prev = head;
        }

        hashLRU() {
            capacity = MAX_INT;
            head = new Node<Key, Value>;
            tail = new Node<Key, Value>;
            head->next = tail;
            tail->prev = head;
        }


        // get the k node in LRU and move it to the head of LRU
        Node<Key, Value> *get(Key k) {
            Node<Key, Value> *node;
            if (map.find(k) != map.end()) {
                node = map[k];
                moveTohead(node);
                return node;
            }// k is existed in map
            else {
                return node;
            }

        }

        // 将node 插入到LRU 的head
        void appendhead(Node<Key, Value> *node) {
            node->prev = head;
            node->next = head->next;
            head->next->prev = node;
            head->next = node;
        }

        // 将 找到的node 移动到head
        void moveTohead(Node<Key, Value> *node) {
            if (node->prev == head) return;
            node->prev->next = node->next;
            node->next->prev = node->prev;
            appendhead(node);
        }

        // put the k node into the head of LRU
        void put(Key k, Value v) {
//            if (k == 17582291756562){
//                cout<< "Jesus, i need You, put 17582291756562" << endl;
//            }
            if (map.find(k) == map.end()) // k is not existed in map
            {
                Node<Key, Value> *node = new Node<Key, Value>(k, v);
                map[k] = node;
                // 判断 size, 如果size = capacity,说明LRU 满了
                if (size == capacity) {
                    poptail();
                }
                size += 1;
                appendhead(node);
            } else {
                map[k]->value = v;
                moveTohead(map[k]);
            }
        }

        //remove the k node from hash LRU
        forceinline void remove(Key k) {
            // 首先找到k 所属的node
//            if (k == 17582291756562){
//                cout<< "Jesus, i need You, remove node 17582291756562" << endl;
//            }
            if (map.find(k) == map.end()) return; // 说明要删除的node 不存在
            Node<Key, Value> *node = map[k];
            node->prev->next = node->next;
            node->next->prev = node->prev;
            map.erase(k);
//            malloc_trim(0);
            size -= 1;
//            delete node;
//            node = NULL;
        }

        //remove the k node from LRU
        forceinline void removenode(Node<Key, Value> *node) {
            // 首先找到k 所属的node
            node->prev->next = node->next;
            node->next->prev = node->prev;
//            if (node->key == 139052311870037896  ){
//                cout<< "Jesus, i need You, delete node 139052311870037896" << endl;
//            }
            if (map.find(node->key) == map.end()){
                cout<< "i need You, my lovely Lord, please come"<< endl;
            }
            else
                map.erase(node->key);
//            malloc_trim(0);
            size -= 1;
//            delete node;
//            node = NULL;
        }

        // pop the tail of the LRU, that the least recent used item
        forceinline Node<Key, Value> *poptail() {
//            map.erase(tail->prev->key);
            map.erase(tail->prev->key);
            Node<Key, Value> *node;
            node = tail->prev;
            tail->prev->prev->next = tail;
            tail->prev = tail->prev->prev;
            size -= 1;
//            if (node->key == 19433408){
//                cout<< "Jesus, i need You, delete 19433408" << endl;
//            }
            return node;
        }

        //get the tail node from local LRU that from this leaf evict key
        Value get_tail() {
            auto tailnode = tail->prev;
//            if (tailnode->key == 19433408){
//                cout<< "Jesus, i need You, delete 19433408" << endl;
//            }
            if (tailnode->value->intrachain.size==1)
                removenode(tailnode);

            return tailnode->value;
        }

        Node<Key, Value>* get_tail1() {
            auto tailnode = tail->prev;
            if (tailnode->value->intrachain.size==1)
                removenode(tailnode);
            return tailnode;
        }

        int deletelru(){

            while (head->next!=NULL){
                auto curnode = head->next;
                curnode->prev = head;
                head->next = curnode->next;
                delete curnode;
            }
            delete head;
//            delete[] tail;
            std::unordered_map<Key, mapvalue>().swap(map);
            malloc_trim(0);
            return 0;
        }

    };

    template<class Key, class Value>
    class localLRU {
    public:
        int size = 0;

        Node<Key, Value> *head;
        Node<Key, Value> *tail;

        localLRU() {
            head = new Node<Key, Value>;
            tail = new Node<Key, Value>;
            head->next = tail;
            tail->prev = head;
        }
        void deletelru(){
            while (head->next!=NULL){
                auto curnode = head->next;
                curnode->prev = head;
                head->next = curnode->next;
//                delete []curnode->value;
//                curnode->value = NULL;
                delete curnode;
                curnode = NULL;
            }
//            malloc_trim(0);
        }

        // 将node 插入到LRU 的head
        forceinline Node<Key, Value> *appendhead(Node<Key, Value> *node) {
            node->prev = head;
            node->next = head->next;
            head->next->prev = node;
            head->next = node;
        }

        // 将 找到的node 移动到head
        forceinline void moveTohead(Node<Key, Value> *node) {
            if (node->prev == head) return;
            node->prev->next = node->next;
            node->next->prev = node->prev;
            appendhead(node);
        }

        // put the k node into the head of LRU
        forceinline Node<Key, Value>*  put(Key k, Value v) {
            Node<Key, Value> *node = new Node<Key, Value>(k, v);
            size += 1;

            appendhead(node);
//            if (k > 4294946590)
//                cout<< "my Lord ,i need You!"<< endl;
            return node;
        }

//        pair<bool,Node<Key, Value>*>  put(Value data) {
//
//            Node<Key, Value> *node = new Node<Key, Value>(data[0], data);
//            size += 1;
//            appendhead(node);
//            return pair<bool,Node<Key, Value>*> (true,node);
//        }

        //remove the k node from local LRU
        forceinline int remove_node(Node<Key, Value> *node) {
            node->prev->next = node->next;
            node->next->prev = node->prev;
//            delete[] node->value;
            delete node;
            node = NULL;
//            malloc_trim(0);
            size -= 1;
            return size;
        }

        // modify the offset (key in lru node), 所有 nodes 中，key > pos 都都需要+ 1
        forceinline void  modify(int pos) {
            // 遍历 intro chain
            // while 循环，直到找到node->key == k;
            Node<Key, Value> *node = head->next;
            while (node != tail)
            {
                if (node->key > pos)
                {
                    node->key += 1;
                }
                node = node->next;
            }
        }



        // pop the tail of the local LRU, that the least recent used item
        forceinline Node<Key, Value> *poptail() {
            Node<Key, Value> *node;
            node = tail->prev;
            tail->prev->prev->next = tail;
            tail->prev = tail->prev->prev;
            size -= 1;
            return node;
        }

    };

    template<class Key>
    class localchain{
    public:
        int size = 0;

        bufNode<Key> *head;
        bufNode<Key> *tail;

        localchain() {
            head = new bufNode<Key>;
            tail = new bufNode<Key>;
            head->next = tail;
            tail->prev = head;
        }
        void deletelru(){
            while (head->next!=NULL){
                auto curnode = head->next;
                curnode->prev = head;
                head->next = curnode->next;
//                delete []curnode->value;
//                curnode->value = NULL;
                delete curnode;
                curnode = NULL;
            }
//            malloc_trim(0);
        }

        // 将node 插入到LRU 的head
        forceinline bufNode<Key> *appendhead(bufNode<Key> *node) {
            node->prev = head;
            node->next = head->next;
            head->next->prev = node;
            head->next = node;
        }

        // 将 找到的node 移动到head
        forceinline void moveTohead(bufNode<Key> *node) {
            if (node->prev == head) return;
            node->prev->next = node->next;
            node->next->prev = node->prev;
            appendhead(node);
        }

        // put the k node into the head of LRU
        forceinline bufNode<Key>*  put(Key k) {
            bufNode<Key> *node = new bufNode<Key>(k);
            size += 1;

            appendhead(node);
//            if (k > 4294946590)
//                cout<< "my Lord ,i need You!"<< endl;
            return node;
        }

//        pair<bool,Node<Key, Value>*>  put(Value data) {
//
//            Node<Key, Value> *node = new Node<Key, Value>(data[0], data);
//            size += 1;
//            appendhead(node);
//            return pair<bool,Node<Key, Value>*> (true,node);
//        }


        //remove the k node from LRU
        forceinline void remove_node(bufNode<Key> *node) {
            node->prev->next = node->next;
            node->next->prev = node->prev;
            delete[] node->value;
            delete node;
//            malloc_trim(0);
            size -= 1;
        }

        // modify the offset (key in lru node), 所有 nodes 中，key > pos 都都需要+ 1
        forceinline void  modify(int pos) {
            // 遍历 intro chain
            // while 循环，直到找到node->key == k;
            bufNode<Key> *node = head->next;
            while (node != tail)
            {
                if (node->key > pos)
                {
                    node->key += 1;
                }
                node = node->next;
            }
        }



        // pop the tail of the local LRU, that the least recent used item
        forceinline bufNode<Key> *poptail() {
            bufNode<Key> *node;
            node = tail->prev;
            tail->prev->prev->next = tail;
            tail->prev = tail->prev->prev;
            size -= 1;
            return node;
        }

    };

}
#endif //EXPERIMENTCC12_FILMADALRU_H
