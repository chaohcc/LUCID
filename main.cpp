
#include <iostream>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>
#include <string.h>
#include <iomanip>
#include <sstream>

#include <sys/time.h>

#include<ctime>

#include "data.h"
#include "workload.h"
#include "film_plus.h"



int test_interleave_insert_query(string filename, double memThreshold, int merge_threshold,
                                 long int datasize, long int numcolumn, int pagesize,unsigned int errbnd,double insert_frac,
                                 unsigned int inner_rebuild_threshold, int argc, char **argv){
    std::cout << "my Lord, i need You, i trust in You! test film interleave insert (append&out-of-order insert)and query" << std::endl;
    struct timeval t1, t2;
    double timeuse;
    double zipf = 0.75;   //zipfian factor
//    unsigned long int numkey = 10000000 ; //ceil(double(datasize*1024*1024)/numcolumn/8);  //int datanum
    unsigned long int numkey = ceil(double(datasize*1024*1024)/numcolumn/8);  //int datanum
    long long test_size = ceil(numkey*2);   // num_key and actual_numkey, the former is the at least number to be inserted to the index. the latter is the keys to guarantee the numkey

    Workload <key_type,key_type*> work;
    work.parse_args(argc, argv);

    // 只读取 文件的前半部分，这种读取形式，如果数据存储偏斜，就会不能很好地用前半部分数据 fit data distribution。
    // 但能够使得 读取文件速度加快
    work.table_size = test_size;

    work.init_table_size = numkey*0.75;
    work.zipf_factor = zipf;

    work.load_keys(test_size);  //  这一步会得到 keys 和 initkeys, keys 是shuffle 过的，initkeys 是sort 过的。

    double reserveMem = 10;
    memThreshold -= reserveMem;

    std::cout << "my Lord, i need You, i trust in You!!!!!" << std::endl;
    gettimeofday(&t1, NULL);

    cout<<"the data set is "<<filename<<",  the data size is "<<datasize<<"M"<< "  numkey is "<<numkey<<endl;
    cout<<"the number of keys is "<<test_size<<", the record size is "<<numcolumn<<endl;
    cout<<"the page size is "<<pagesize*8<<",  the available memory is "<< memThreshold<< endl;

    unsigned long init_num_key = ceil(numkey*0.75);




    filminsert::test_interleave_insert_query(errbnd,numkey,pagesize,filename,memThreshold,
                                             reserveMem,numcolumn, merge_threshold,work,test_size,datasize,
                                             zipf,init_num_key, inner_rebuild_threshold);

    delete work.keys;
    std::vector<key_type>().swap(work.init_keys);

    gettimeofday(&t2, NULL);
    timeuse = (t2.tv_sec - t1.tv_sec) + (double) (t2.tv_usec - t1.tv_usec) / 1000000.0;

    cout << "able o_direct disk access time = " << timeuse << endl;  //输出时间（单位：ｓ）
    return 0;
}



int test_inner_rebuild_threshold(string filename, double memThreshold, long int datasize, long int numcolumn, int pagesize,unsigned short int errbnd, int merge_threshold,
                                 double insert_ratio,int argc, char **argv){
    unsigned int inner_rebuild_thresholds[] = {16,32,64,128};
    for (int k = 0; k < end(inner_rebuild_thresholds)-begin(inner_rebuild_thresholds);k++){
        unsigned int inner_rebuild_threshold = inner_rebuild_thresholds[k] ;
        test_interleave_insert_query(filename, memThreshold, merge_threshold,
                                     datasize, numcolumn, pagesize,errbnd,insert_ratio,inner_rebuild_threshold,argc, argv);

    }
    return 0;
}

int test_merge_threshold(string filename, double memThreshold, long int datasize, long int numcolumn, int pagesize,unsigned short int errbnd, double insert_ratio,
                         unsigned int inner_rebuild_threshold,int argc, char **argv){
    double merge_ratios[] = {0.5} ;  // 0.1,0.2,0.3,0.4,
    for (int i = 0; i < (end(merge_ratios)-begin(merge_ratios));i++) {

        double ratio = merge_ratios[i];
        int merge_threshold = ceil(ratio * pagesize/numcolumn);
        test_interleave_insert_query(filename,memThreshold, merge_threshold, datasize, numcolumn, pagesize,errbnd, insert_ratio,inner_rebuild_threshold,argc, argv);
    }
    return 0;
}


int test_record(string filename, double memThreshold, long int datasize, int merge_threshold, int pagesize,unsigned short int errbnd, double insert_ratio,
                unsigned int inner_rebuild_threshold,int argc, char **argv){
    int merge_ratios[] = {128};  // 8，16，32，64
    for (int i = 0; i < (end(merge_ratios)-begin(merge_ratios));i++) {

        int numcolumn = merge_ratios[i];
        merge_threshold = ceil(0.5*pagesize/numcolumn);
        test_interleave_insert_query(filename,memThreshold, merge_threshold, datasize, numcolumn, pagesize,errbnd, insert_ratio,inner_rebuild_threshold,argc, argv);
    }
    return 0;
}



int test_mem(string filename, int numcolumn, long int datasize, int merge_threshold, int pagesize,unsigned short int errbnd, double insert_ratio,
             unsigned int inner_rebuild_threshold,int argc, char **argv){
    double merge_ratios[] = {2048};  // 0.1,0.2,0.3,0.4,
    for (int i = 0; i < (end(merge_ratios)-begin(merge_ratios));i++) {

        double memThreshold = merge_ratios[i];

        test_interleave_insert_query(filename,memThreshold, merge_threshold, datasize, numcolumn, pagesize,errbnd, insert_ratio,inner_rebuild_threshold,argc, argv);
    }
    return 0;
}


int test_error(string filename, double memThreshold, long int datasize, long int numcolumn, int pagesize,int merge_threshold, double insert_ratio,
               unsigned int inner_rebuild_threshold,int argc, char **argv){
    int merge_ratios[] = {4,8,16,32};  // 4,8,16,32   64,128,256
    for (int i = 0; i < (end(merge_ratios)-begin(merge_ratios));i++) {

        int errbnd = merge_ratios[i];

        test_interleave_insert_query(filename,memThreshold, merge_threshold, datasize, numcolumn, pagesize,errbnd, insert_ratio,inner_rebuild_threshold,argc, argv);
    }
    return 0;
}


int test_inner_threshold(string filename, double memThreshold, long int datasize, long int numcolumn, int pagesize,unsigned short int errbnd, double insert_ratio,
                         int merge_threshold,int argc, char **argv){
    double inner_thresholds[] = {32};  // 8,32,128,512,2048,8192,16384
    for (int i = 0; i < (end(inner_thresholds)-begin(inner_thresholds));i++) {


        int inner_rebuild_threshold = inner_thresholds[i];
        test_interleave_insert_query(filename,memThreshold, merge_threshold, datasize, numcolumn, pagesize,errbnd, insert_ratio,inner_rebuild_threshold,argc, argv);
    }
    return 0;
}


int test_pagesize(string filename, double memThreshold, long int datasize, long int numcolumn, unsigned short int errbnd, double insert_ratio,
                  unsigned int inner_rebuild_threshold,int argc, char **argv,int merge_threshold){
    int pagesizes[] = {4,8,16,32,64,128};
    for (int i = 0; i < (end(pagesizes)-begin(pagesizes));i++) {

        int pagesize = pagesizes[i]*1024/8;
        int merge_threshold = ceil(0.5 * pagesize/numcolumn);
        test_interleave_insert_query(filename,memThreshold, merge_threshold, datasize, numcolumn, pagesize,errbnd, insert_ratio,inner_rebuild_threshold,argc, argv);
    }
    return 0;
}

unsigned int choose_inner_rebuild(string workload){
    unsigned int inner_rebuild_thresholds[] = {256,512,1024,1536,2048};
    unsigned int inner_rebuild_threshold = 0;
    if (workload == "zipf")
        inner_rebuild_threshold = inner_rebuild_thresholds[0];
    else if (workload == "zipfrandom")
        inner_rebuild_threshold = inner_rebuild_thresholds[1];
    else if (workload == "hotspot")
        inner_rebuild_threshold = inner_rebuild_thresholds[2];
    else
        inner_rebuild_threshold = inner_rebuild_thresholds[2];

    return inner_rebuild_threshold;
}

int main() {

    std::cout << "my Lord, i need You, i trust in You! start the experiment " << std::endl;

    unsigned int errbnd = 16;
    double memThreshold = 8192;
    int pagesize = 1024*64/8;
    string filename = "wiki_ts";  //queryname
    long int numcolumn = 16;
    long int datasize = 4096;
    double zipf = 0.5;
    string workload = "zipf";
    datasize = 2048*2;
    memThreshold = 1024*2;

    double insert_ratio = 0.5;
    int merge_threshold = ceil(0.5*pagesize/numcolumn);
    unsigned int inner_rebuild_threshold = 256;
    errbnd = 16 ;

    // out-of-order insertion 的代码片段
    datasize = 1024*4;
    memThreshold = 1024*2; // available memory
    filename = "wise";
//    std::string datapath = "--keys_file=/home/wamdm/chaohong/TLIdata/keydata/wise";
    std::string datapath = "--keys_file=/data/TLIdata/keydata/wise";

    char *datafile = const_cast<char *>(datapath.c_str());
    if (filename == "fb")
        errbnd = 64;
    else if (filename == "planet")
        errbnd = 32;
    else
        errbnd = 16;

//    sample_distribution = zipf, uniform, hotspot, zipfrandom

//    --table_size= -1 表示读取全部的数据
//    "--update=0.1" ,"--scan=0.1" ,  "--insert=0.4"
//    --keys_file=/home/wamdm/chaohong/clionDir/ALEX_1D_binary/wise
//    /data/chaohong/ALEX_1D_binary/genome
//   binary longitudes binary double  longitudes-200M.bin.data
//    binary  unsigned long osm_cellids_800M_uint64.zst, osm
// binary fb unsigned long  fb_200M_uint64.zst
// binary wiki_ts unsigned long wiki_ts_200M_uint64
// binary books unsigned long  books_200M_uint32
// binary stack unsigned long  stack

//"--scan=0.1",
//    "--buffer=0.4",
//    "--insert=0.5","--read=0.1","--scan=0.4"

    pagesize = 1024*32/8;
//    "--insert=0.5","--read=0.4","--scan=0.1"
    merge_threshold = ceil(0.5*pagesize/numcolumn);

    // big data size
//    char* argv1_1[] = {"test","--keys_file=/data/chaohong/ALEX_1D_binary/ycsb_bigsize", "--keys_file_type=text","--batch_num=2","--scan_num=0",
//                     "--insert=0","--read=1","--scan=0","--sample_distribution=zipf","--buffer=0",
//                     "--operations_num=1000000","--table_size=-1","--init_table_ratio=1" };
//    int argc1_1 = (end(argv1_1)-begin(argv1_1));
//    test_merge_threshold("ycsb_bigsize", 64*1024, 64*1024, numcolumn, pagesize,errbnd,insert_ratio,inner_rebuild_threshold,argc1_1, argv1_1) ;


    char* argv3[] = {"test",datafile, "--keys_file_type=binary","--batch_num=21","--scan_num=0",
                     "--insert=0.4","--read=0.4","--scan=0","--delete=0.2","--sample_distribution=random","--buffer=0.3",
                     "--operations_num=300000","--table_size=-1","--init_table_ratio=0.75" };


    int argc3 = (end(argv3)-begin(argv3));
//    test_error("osm", memThreshold, datasize, numcolumn, pagesize,merge_threshold,insert_ratio,inner_rebuild_threshold,argc3, argv3);
//    test_mem(filename, numcolumn, datasize, merge_threshold, pagesize,errbnd,insert_ratio,inner_rebuild_threshold,argc3, argv3) ;
//    test_pagesize(filename,  memThreshold, datasize, numcolumn, errbnd,insert_ratio,inner_rebuild_threshold,argc3, argv3,merge_threshold) ;
    test_merge_threshold(filename, memThreshold, datasize, numcolumn, pagesize,errbnd,insert_ratio,inner_rebuild_threshold,argc3, argv3) ;
//    test_merge_threshold(filename,  memThreshold, datasize, numcolumn, pagesize,errbnd,insert_ratio,inner_rebuild_threshold,argc3, argv3) ;
//    test_record(filename, memThreshold, datasize, merge_threshold, pagesize,errbnd,insert_ratio,inner_rebuild_threshold,argc3, argv3) ;
//    test_inner_threshold(filename,  memThreshold, datasize, numcolumn, pagesize,errbnd,insert_ratio,merge_threshold,argc3, argv3);

//    test_inner_threshold(filename,  memThreshold, datasize, numcolumn, pagesize,errbnd,insert_ratio,merge_threshold,argc1, argv1);



    char* argv12[] = {"test",datafile, "--keys_file_type=binary","--batch_num=21","--scan_num=0",
                      "--insert=0.64","--read=0.16","--scan=0","--delete=0.2","--sample_distribution=random","--buffer=0.3",
                      "--operations_num=300000","--table_size=-1","--init_table_ratio=0.75" };
    pagesize = 1024*32/8;
    int argc12 = (end(argv12)-begin(argv12));
    test_merge_threshold(filename, memThreshold, datasize, numcolumn, pagesize,errbnd,insert_ratio,inner_rebuild_threshold,argc12, argv12) ;
//    test_pagesize(filename,  memThreshold, datasize, numcolumn, 16,insert_ratio,inner_rebuild_threshold,argc12, argv12,merge_threshold) ;
//    test_error("fb", memThreshold, datasize, numcolumn, pagesize,merge_threshold,insert_ratio,inner_rebuild_threshold,argc12, argv12) ;

    char* argv121[] = {"test",datafile, "--keys_file_type=binary","--batch_num=21","--scan_num=0",
                      "--insert=0.16","--read=0.64","--scan=0","--delete=0.2","--sample_distribution=random","--buffer=0.3",
                      "--operations_num=300000","--table_size=-1","--init_table_ratio=0.75" };
    pagesize = 1024*32/8;
    int argc121 = (end(argv121)-begin(argv121));
    test_merge_threshold(filename, memThreshold, datasize, numcolumn, pagesize,errbnd,insert_ratio,inner_rebuild_threshold,argc121, argv121) ;

    char* argv1[] = {"test",datafile, "--keys_file_type=binary","--batch_num=21","--scan_num=0",
                     "--insert=0.12","--read=0.48","--scan=0","--delete=0.4","--sample_distribution=zipf","--buffer=0.3",
                     "--operations_num=300000","--table_size=-1","--init_table_ratio=0.75" };
    pagesize = 1024*32/8;
    int argc1 = (end(argv1)-begin(argv1));
//    test_record(filename, memThreshold, datasize, merge_threshold, pagesize,errbnd,insert_ratio,inner_rebuild_threshold,argc1, argv1) ;
//    test_error("wise", memThreshold, datasize, numcolumn, pagesize,merge_threshold,insert_ratio,inner_rebuild_threshold,argc1, argv1) ;
//    test_merge_threshold(filename, memThreshold, datasize, numcolumn, pagesize,errbnd,insert_ratio,inner_rebuild_threshold,argc1, argv1) ;
//    test_pagesize(filename,  memThreshold, datasize, 16, errbnd,insert_ratio,inner_rebuild_threshold,argc1, argv1,merge_threshold) ;



// inner rebuild
    char* argv11[] =  {"test",datafile, "--keys_file_type=binary","--batch_num=21","--scan_num=0",
                       "--insert=0.3","--read=0.3","--scan=0","--delete=0.4","--sample_distribution=zipf","--buffer=0.3",
                       "--operations_num=300000","--table_size=-1","--init_table_ratio=0.75" };

    pagesize = 1024*32/8;
    int argc11 = (end(argv11)-begin(argv11));
//    test_merge_threshold2(filename, memThreshold, datasize, numcolumn, pagesize,errbnd,insert_ratio,inner_rebuild_threshold,argc11, argv11) ;
    //test_merge_threshold(filename, memThreshold, datasize, numcolumn, pagesize,errbnd,insert_ratio,inner_rebuild_threshold,argc11, argv11) ;
//    test_inner_threshold(filename,  memThreshold, datasize, numcolumn, pagesize,errbnd,insert_ratio,merge_threshold,argc11, argv11);
//    test_record(filename, memThreshold, datasize, merge_threshold, pagesize,errbnd,insert_ratio,inner_rebuild_threshold,argc11, argv11) ;
//    test_mem(filename, numcolumn, datasize, merge_threshold, pagesize,errbnd,insert_ratio,inner_rebuild_threshold,argc11, argv11) ;
//    test_error("planet", memThreshold, datasize, numcolumn, pagesize,merge_threshold,insert_ratio,inner_rebuild_threshold,argc11, argv11) ;
//    test_pagesize(filename,  memThreshold, datasize, numcolumn, errbnd,insert_ratio,inner_rebuild_threshold,argc11, argv11,merge_threshold) ;


////    "--insert=0.5","--read=0.25","--scan=0.25",

// test range size   "--operations_num=500000"
//    char* argv2[] = {"test","--keys_file=/data/chaohong/ALEX_1D_binary/genome", "--keys_file_type=binary","--batch_num=3","--scan_num=100000",
//                     "--insert=0.05","--read=0","--scan=0.95","--sample_distribution=zipf","--buffer=0.3",
//                     "--operations_num=500000","--table_size=-1","--init_table_ratio=0.75" };


    char* argv2[] = {"test","--keys_file=/data/chaohong/ALEX_1D_binary/covid", "--keys_file_type=binary","--batch_num=31","--scan_num=100",
                     "--insert=0.5","--read=0.5","--scan=0","--sample_distribution=zipf","--buffer=0.15",
                     "--operations_num=300000","--table_size=-1","--init_table_ratio=0.75" };
    pagesize = 1024*32/8;
    int argc2 = (end(argv2)-begin(argv2));
//    test_error("covid", memThreshold, datasize, numcolumn, pagesize,merge_threshold,insert_ratio,inner_rebuild_threshold,argc2, argv2) ;
//    test_merge_threshold("genome", memThreshold, datasize, numcolumn, pagesize,errbnd,insert_ratio,inner_rebuild_threshold,argc2, argv2) ;
//      test_pagesize(filename,  memThreshold, datasize, numcolumn, errbnd,insert_ratio,inner_rebuild_threshold,argc2, argv2,merge_threshold) ;
//    test_merge_threshold(filename, memThreshold, datasize, numcolumn, pagesize,errbnd,insert_ratio,inner_rebuild_threshold,argc2, argv2) ;
////    test_merge_threshold(filename, memThreshold, datasize, numcolumn, pagesize,errbnd,insert_ratio,64,argc2, argv2) ;
////    test_merge_threshold(filename, memThreshold, datasize, numcolumn, pagesize,errbnd,insert_ratio,256,argc2, argv2) ;
////    test_merge_threshold(filename, memThreshold, datasize, numcolumn, pagesize,errbnd,insert_ratio,1024,argc2, argv2) ;

    pagesize = 1024*32/8;
    merge_threshold = ceil(0.5*pagesize/numcolumn);
    char* argv13[] = {"test","--keys_file=/data/chaohong/ALEX_1D_binary/lognormal", "--keys_file_type=binary","--batch_num=31","--scan_num=100",
                      "--insert=0.5","--read=0.5","--scan=0","--sample_distribution=zipf","--buffer=0.15",
                      "--operations_num=300000","--table_size=-1","--init_table_ratio=0.75" };

    int argc13 = (end(argv13)-begin(argv13));
//    test_pagesize(filename,  memThreshold, datasize, numcolumn, 64,insert_ratio,inner_rebuild_threshold,argc13, argv13,merge_threshold) ;
//    test_merge_threshold("longitudes",  memThreshold, datasize, numcolumn, pagesize,errbnd,insert_ratio,inner_rebuild_threshold,argc13, argv13) ;
//    test_error("lognormal", memThreshold, datasize, numcolumn, pagesize,merge_threshold,insert_ratio,inner_rebuild_threshold,argc13, argv13) ;

    char* argv103[] = {"test",datafile, "--keys_file_type=binary","--batch_num=31","--scan_num=100",
                       "--insert=0.2","--read=0.8","--scan=0","--sample_distribution=zipf","--buffer=0.2",
                       "--operations_num=300000","--table_size=-1","--init_table_ratio=0.75" };

    int argc103 = (end(argv103)-begin(argv103));
//    test_pagesize(filename,  memThreshold, datasize, numcolumn, 64,insert_ratio,inner_rebuild_threshold,argc13, argv13,merge_threshold) ;
//    test_merge_threshold(filename,  memThreshold, datasize, numcolumn, pagesize,errbnd,insert_ratio,inner_rebuild_threshold,argc103, argv103) ;
//
    char* argv113[] = {"test",datafile, "--keys_file_type=binary","--batch_num=31","--scan_num=100",
                       "--insert=0.2","--read=0.8","--scan=0","--sample_distribution=zipf","--buffer=0.2",
                       "--operations_num=300000","--table_size=-1","--init_table_ratio=0.75" };

    int argc113 = (end(argv13)-begin(argv13));
//    test_pagesize(filename,  memThreshold, datasize, numcolumn, 64,insert_ratio,inner_rebuild_threshold,argc13, argv13,merge_threshold) ;
//    test_merge_threshold(filename,  memThreshold, datasize, numcolumn, pagesize,errbnd,insert_ratio,inner_rebuild_threshold,argc113, argv113) ;
//

    char* argv123[] = {"test",datafile, "--keys_file_type=binary","--batch_num=31","--scan_num=100",
                       "--insert=0.2","--read=0.8","--scan=0","--sample_distribution=zipf","--buffer=0.2",
                       "--operations_num=300000","--table_size=-1","--init_table_ratio=0.75" };

    int argc123 = (end(argv13)-begin(argv13));
//    test_pagesize(filename,  memThreshold, datasize, numcolumn, 64,insert_ratio,inner_rebuild_threshold,argc13, argv13,merge_threshold) ;
//    test_merge_threshold(filename,  memThreshold, datasize, numcolumn, pagesize,errbnd,insert_ratio,inner_rebuild_threshold,argc123, argv123) ;
//

    char* argv133[] = {"test",datafile, "--keys_file_type=binary","--batch_num=31","--scan_num=100",
                       "--insert=0.2","--read=0.8","--scan=0","--sample_distribution=zipf","--buffer=0.4",
                       "--operations_num=300000","--table_size=-1","--init_table_ratio=0.75" };

    int argc133 = (end(argv13)-begin(argv13));
//    test_pagesize(filename,  memThreshold, datasize, numcolumn, 64,insert_ratio,inner_rebuild_threshold,argc13, argv13,merge_threshold) ;
//    test_merge_threshold(filename,  memThreshold, datasize, numcolumn, pagesize,errbnd,insert_ratio,inner_rebuild_threshold,argc133, argv133) ;
//

    char* argv1331[] = {"test",datafile, "--keys_file_type=binary","--batch_num=31","--scan_num=100",
                        "--insert=0.5","--read=0.5","--scan=0","--sample_distribution=zipf","--buffer=0.1",
                        "--operations_num=300000","--table_size=-1","--init_table_ratio=0.75" };

    int argc1331 = (end(argv13)-begin(argv13));
//    test_pagesize(filename,  memThreshold, datasize, numcolumn, 64,insert_ratio,inner_rebuild_threshold,argc13, argv13,merge_threshold) ;
//    test_merge_threshold(filename,  memThreshold, datasize, numcolumn, pagesize,errbnd,insert_ratio,inner_rebuild_threshold,argc1331, argv1331) ;

    char* argv1332[] = {"test",datafile, "--keys_file_type=binary","--batch_num=31","--scan_num=100",
                        "--insert=0.5","--read=0.5","--scan=0","--sample_distribution=zipf","--buffer=0.2",
                        "--operations_num=300000","--table_size=-1","--init_table_ratio=0.75" };

    int argc1332 = (end(argv13)-begin(argv13));
//    test_pagesize(filename,  memThreshold, datasize, numcolumn, 64,insert_ratio,inner_rebuild_threshold,argc13, argv13,merge_threshold) ;
//    test_merge_threshold(filename,  memThreshold, datasize, numcolumn, pagesize,errbnd,insert_ratio,inner_rebuild_threshold,argc1332, argv1332);

    char* argv1333[] = {"test",datafile, "--keys_file_type=binary","--batch_num=31","--scan_num=100",
                        "--insert=0.5","--read=0.5","--scan=0","--sample_distribution=zipf","--buffer=0.3",
                        "--operations_num=300000","--table_size=-1","--init_table_ratio=0.75" };

    int argc1333 = (end(argv13)-begin(argv13));
//    test_pagesize(filename,  memThreshold, datasize, numcolumn, 64,insert_ratio,inner_rebuild_threshold,argc13, argv13,merge_threshold) ;
//    test_merge_threshold(filename,  memThreshold, datasize, numcolumn, pagesize,errbnd,insert_ratio,inner_rebuild_threshold,argc1333, argv1333);

    char* argv1334[] = {"test",datafile, "--keys_file_type=binary","--batch_num=31","--scan_num=100",
                        "--insert=0.8","--read=0.2","--scan=0","--sample_distribution=zipf","--buffer=0.1",
                        "--operations_num=300000","--table_size=-1","--init_table_ratio=0.75" };

    int argc1334 = (end(argv13)-begin(argv13));
//    test_pagesize(filename,  memThreshold, datasize, numcolumn, 64,insert_ratio,inner_rebuild_threshold,argc13, argv13,merge_threshold) ;
//    test_merge_threshold(filename,  memThreshold, datasize, numcolumn, pagesize,errbnd,insert_ratio,inner_rebuild_threshold,argc1334, argv1334);

    return 0;
}


/*  // 对于int 类型的值，分高低字节编辑和读取
unsigned int iTest=0;
unsigned short int *piTest=(unsigned short int *)&iTest;
*piTest=1024;	//低16位值
piTest++;
*piTest=768;
vector<unsigned int> arr = {100,6778,12345};
for (int i = 0; i < 3; i ++){
    auto dual_off_iter = arr.end()-1 ;
    cout <<"the high 16 bits " << ((*dual_off_iter)>>16) << " the low 16 bits "<< ((*dual_off_iter)&0xFFFF)<< endl;
    unsigned short int *piTest=(unsigned short int *)&(*dual_off_iter);
    *piTest=i+1;	//低16位值
    cout <<"the high 16 bits " << ((*dual_off_iter)>>16) << " the low 16 bits "<< ((*dual_off_iter)&0xFFFF)<< endl;
    cout << "i need you, my Lord! thank You!" << endl;
}
*/


