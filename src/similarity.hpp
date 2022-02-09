#ifndef SIMILARITYINDEX
#define SIMILARITYINDEX

#include "internal.hpp"
#include <iostream>
#include <map>
#include <omp.h>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <time.h>

using namespace std;

//パラメーターセッティング
#define CONSTANT_FOR_RANK_CALC 10 // Rank計算用定数
#define CSD_SET_CRITERIA 1        // VSIDSの値がx未満であれば価値がないと判断する。0=VSIDSに何かしらの値が入っていればOK。1=直近の学習説に含まれているものにほぼ限定される
#define RESTART_LIMIT 100000
#define CHANGE_SCORE_RATIO 0.1
#define CHANGE_SCORE_INCRE 10000
#define ALPHA_TO_JUDGE_SSI 2
#define LIMIT_SAVING_SSI 1000
#define LIMIT_SAVING_CSD 1000
#define CHANGE_INTERVAL 10
#define STABLE_ONLY_MODE false

struct csd_element
{
    int rank;
    bool phase;
    double value;
};

struct CSD
{
    int nonZeroVars;
    vector<csd_element> data;
};

enum similarityLevel
{
    high,
    normal,
    low
};

extern CSD get_CSD(vector<double>, vector<int>, bool, CaDiCaL::Phases);
extern double calculate_SSI(CSD, CSD);
extern void save_CSD(CSD);
extern CSD get_prevCSD(int);
extern CSD tmp_csd, saved_csd;
extern vector<int> set_qtab(CaDiCaL::Queue, CaDiCaL::Links);

extern vector<double> SSI_database;
extern bool check_ssi_table(int thread_num);
extern vector<double> change_search_space(vector<double> &score_table, CaDiCaL::ScoreSchedule &scores, double scoreInc);
similarityLevel judge_SSI_score(double ssi);

CSD init_csd(size_t var_size);

#endif
