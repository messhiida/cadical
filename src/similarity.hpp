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
#define ALPHA_TO_JUDGE_SSI 2
#define CHANGE_SCORE_RATIO 0.1
#define CHANGE_SCORE_INCRE 10000
#define ACTION_THREASHOLD 2
//#define RESTART_LIMIT 10000
#define CONSTANT_FOR_RANK_CALC 10
#define CSD_SCORE_CRITERIA 0 //VSIDSの値がx未満であれば価値がないと判断する。0=VSIDSに何かしらの値が入っていればOK。1=直近の学習説に含まれているものにほぼ限定される
#define CSD_VALID_VAR_RATIO 0.333
#define LIMIT_SAVING_SSI 1000
#define LIMIT_SAVING_CSD 101

enum similarityLevel
{
    high,
    normal,
    low
};

extern vector<double> SSI_database;
extern vector<vector<array<double, 3>>> csd_database;
extern bool check_ssi_table(int thread_num);
extern vector<double> change_search_space(vector<double> &score_table, CaDiCaL::ScoreSchedule &scores, double scoreInc);
vector<array<double, 3>> get_CSD(vector<double> scoreTable, vector<signed char> phases, CaDiCaL::ScoreSchedule scores);
double calculate_SSI(vector<array<double, 3>> csd1, vector<array<double, 3>> csd2);
void save_CSD(vector<array<double, 3>> csd);
similarityLevel judge_SSI_score(double ssi);

#endif
