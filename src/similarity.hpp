#ifndef SIMILARITYINDEX
#define SIMILARITYINDEX

#include <iostream>
#include <map>
#include "internal.hpp"

using namespace std;

#define RESTART_LIMIT 5000
#define CONSTANT_FOR_RANK_CALC 10
#define ALPHA_TO_JUDGE_SSI 1

#define LIMIT_SAVING_SSI 100
extern vector<double> SSI_database;
#define LIMIT_SAVING_CSD 101
extern vector< map<int, vector<double> > > csd_database;

extern int scoreCriteria; 
extern double affectRatioForChangingSearch; //下位何％をpopするか
extern double incrementalNumForChangingSearch; //10000回分のVSIDS bumpに相当する分の変更を加える
enum similarityLevel {high, normal, low};

void showResult(vector<double> s);
void draw_VSIDS_graph (vector<double> s);

double countValidScoreVariables(vector<double> scores);

map<int, vector<double> > get_CSD (vector<double> scores, vector<signed char> phases);
double calculate_SSI (map<int, vector<double> > csd1, map<int, vector<double> > csd2);

void save_SSI (double ssi);
void save_CSD (map<int, vector<double> > csd);
similarityLevel judge_SSI_score(double ssi);
vector<int> read_learntClause(CaDiCaL::Clause*);

//Todo ここまだこれから
vector<double> change_search_space(vector<double> scores, double scoreInc);

#endif