#ifndef SIMILARITYINDEX
#define SIMILARITYINDEX

#include <iostream>
#include <map>
#include "internal.hpp"
#include <omp.h>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>

using namespace std;

//パラメーターセッティング
#define RESTART_LIMIT 10000
#define CONSTANT_FOR_RANK_CALC 10
#define ALPHA_TO_JUDGE_SSI 1
#define CSD_SCORE_CRITERIA 1 //VSIDSの値がx未満であれば価値がないと判断する。0=VSIDSに何かしらの値が入っていればOK。1=直近の学習説に含まれているものにほぼ限定される
#define LIMIT_SAVING_SSI 100
#define LIMIT_SAVING_CSD 101
#define FIXED_DECISION_AT_TOP true
#define PARALLEL_NUM 3
#define PARA_SHARED_CLAUSE_LBD 2 //これ以下は共有する

enum similarityLevel
{
    high,
    normal,
    low
};

extern vector<double> SSI_database;
extern vector<map<int, vector<double>>> csd_database;
extern int sharedData;
extern bool para_finished;
extern int set_parallel_seed(int thread_num);
extern void submit_shared_learntClause(int thread_num, CaDiCaL::Clause *lc);
extern vector<CaDiCaL::Clause *> import_shared_learntClause();

map<int, vector<double>> get_CSD(vector<double> scores, vector<signed char> phases);
double calculate_SSI(map<int, vector<double>> csd1, map<int, vector<double>> csd2);
void save_CSD(map<int, vector<double>> csd);

similarityLevel judge_SSI_score(double ssi);
vector<int> convert_learntClause_to_vector(CaDiCaL::Clause *);

extern vector<int> getSameLearntClause_restart(string &outputFile);

//Todo ここまだこれから
vector<double> change_search_space(vector<double> scores, double scoreInc);

#endif