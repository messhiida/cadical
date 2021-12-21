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
#define PARALLEL_NUM 2
#define ALPHA_TO_JUDGE_SSI 3
#define CHANGE_SCORE_RATIO 0.1
#define CHANGE_SCORE_INCRE 10000
#define ACTION_THREASHOLD 2
//#define RESTART_LIMIT 10000
//#define FIXED_DECISION_AT_TOP true
#define PARA_SHARED_CLAUSE_LBD 3 //これ以下は共有する
#define CONSTANT_FOR_RANK_CALC 10
#define CSD_SCORE_CRITERIA 0 //VSIDSの値がx未満であれば価値がないと判断する。0=VSIDSに何かしらの値が入っていればOK。1=直近の学習説に含まれているものにほぼ限定される
#define CSD_VALID_VAR_RATIO 0.333
#define LIMIT_SAVING_SSI 1000
#define LIMIT_SAVING_CSD 101
#define LIMIT_SHARED_CSD 10

enum similarityLevel
{
    high,
    normal,
    low
};

extern vector<double> SSI_database;
extern vector<vector<array<double, 3>>> csd_database;
extern bool para_finished;
extern int set_parallel_seed(int thread_num);
extern void submit_shared_learntClause(int thread_num, vector<CaDiCaL::Clause *> &lc);
extern vector<CaDiCaL::Clause *> import_shared_learntClause(int thread_num, vector<CaDiCaL::Clause *> &clauses, CaDiCaL::Stats &stats);

extern void announce_para_finished();
extern bool check_para_finished();

extern array<vector<array<double, 3>>, PARALLEL_NUM> shared_csd;
extern bool check_ssi_table(int thread_num);
extern void set_bool_to_action_table(int i, int j, bool b);

extern void submit_csd(int thread_num, vector<array<double, 3>> csd);
extern int update_worker_action_table();
extern vector<double> change_search_space(vector<double> &score_table, CaDiCaL::ScoreSchedule &scores, double scoreInc);

vector<array<double, 3>> get_CSD(vector<double> scoreTable, vector<signed char> phases, CaDiCaL::ScoreSchedule scores);
double calculate_SSI(vector<array<double, 3>> csd1, vector<array<double, 3>> csd2);
void save_CSD(vector<array<double, 3>> csd);

similarityLevel judge_SSI_score(double ssi);
vector<int> convert_learntClause_to_vector(CaDiCaL::Clause *);

extern vector<int> getSameLearntClause_restart(string &outputFile);

#endif
