#include "internal.hpp"

using namespace std;
vector<double> SSI_database;
bool para_finished = false;
array<vector<array<double, 3>>, PARALLEL_NUM> shared_csd;
array<vector<CaDiCaL::Clause *>, PARALLEL_NUM> shared_learntClause;

bool check_ssi_table(int thread_num)
{
    for (int i = thread_num; i < PARALLEL_NUM; i++)
    {
        if (i == thread_num)
            continue;

        vector<array<double, 3>> my_csd, comp_csd;
#pragma omp critical(SHARED_CSD)
        {
            my_csd = shared_csd[thread_num];
            comp_csd = shared_csd[i];
        }
        double ssi = calculate_SSI(my_csd, comp_csd);
        similarityLevel sl = judge_SSI_score(ssi);
        if (sl == high)
        {
            printf("Similarity high: [%d][%d]\n", thread_num, i);
            return true;
        }
    }
    return false;
}

void submit_csd(int thread_num, vector<array<double, 3>> csd)
{
#pragma omp critical(SHARED_CSD)
    shared_csd[thread_num] = csd;
}

void announce_para_finished()
{
#pragma omp critical(PARA_FINISHED)
    para_finished = true;
}

bool check_para_finished()
{
    bool tmp_fin;
#pragma omp critical(PARA_FINISHED)
    tmp_fin = para_finished;
    if (tmp_fin == true)
        return true;
    else
        return false;
}

vector<double> change_search_space(vector<double> &score_table, CaDiCaL::ScoreSchedule &scores, double inc)
{
    double incrementalScore = CHANGE_SCORE_INCRE * inc; //score_incのINCRE回分をup (CIRと同じであれば10000回）
    while (incrementalScore > 1e150)                    //double上限には引っかからないとは思うが、念の為
        incrementalScore /= 10;

    vector<int> tmp_lit;
    auto it = scores.end();
    double limit = score_table.size() * CHANGE_SCORE_RATIO;
    double counter = 0;
    //一度scoreを元に回して置かなければ、updateにおいてscoreのiteratorが崩れてしまうため2回ループを回す
    while (it != scores.begin())
    {
        --it; //一番はじめにデクリメントしてやる必要あり
        if (counter >= limit)
            break; //一定数以上の変数を触るとBreak
        tmp_lit.push_back(*it);
        counter++;
    }
    for (auto l : tmp_lit)
    {
        score_table[l] += incrementalScore;
        scores.update(l);
    }

    return score_table;
}

void submit_shared_learntClause(int thread_num, vector<CaDiCaL::Clause *> &lc)
{
    for (CaDiCaL::Clause *c : lc)
    {
        for (int i = 0; i < PARALLEL_NUM; i++)
        {
            if (i != thread_num) //thread_num == 自分自身 の為skip
            {
#pragma omp critical(SHARED_CLAUSE)
                shared_learntClause[i].push_back(c);
            }
        }
    }
    lc.clear();
}

vector<CaDiCaL::Clause *> import_shared_learntClause(int thread_num, vector<CaDiCaL::Clause *> &clauses, CaDiCaL::Stats &stats)
{
    vector<CaDiCaL::Clause *> slc; //shared learnt clause

#pragma omp critical(SHARED_CLAUSE)
    slc = shared_learntClause[thread_num];

    for (const auto &c : slc)
    {
        int size = c->size;
        stats.learned.literals += size;
        stats.learned.clauses++;

        stats.units += (size == 1);
        stats.binaries += (size == 2);

        stats.added.total++;
        stats.current.total++;
        stats.added.total++;

        if (c->redundant)
        {
            stats.current.redundant++;
            stats.added.redundant++;
        }
        else
        {
            stats.irrbytes += c->bytes();
            stats.current.irredundant++;
            stats.added.irredundant++;
        }
        //TODO ここと
        clauses.push_back(c);
        stats.learned.literals += c->size;
        stats.learned.clauses++;
    }

#pragma omp critical(SHARED_CLAUSE)
    shared_learntClause[thread_num].clear();

    return slc;
}

int set_parallel_seed(int thread_num)
{
    int para_seed = 0;
    para_seed = 2048 * thread_num; //2048はmagic number
    return para_seed;
}

double count_validScoreVars(vector<double> stab)
{
    double count = 1;
    for (double s : stab)
        if (s > CSD_SCORE_CRITERIA)
            count++;
    return count;
}

vector<array<double, 3>> get_CSD(vector<double> scoreTable, vector<signed char> phases, CaDiCaL::ScoreSchedule scores)
{
    int var_size = (int)scoreTable.size();
    vector<array<double, 3>> csd(var_size); //csd[var] = {rank, phase, value},　不足分=下で定義されない分はzeroで埋められる
    double rank = 0.0;

    for (auto it = scores.begin(); it != scores.end(); ++it)
    {
        rank++;

        unsigned var_index = *it;
        double score = scoreTable[*it];
        if (score <= (double)CSD_SCORE_CRITERIA)
            break;

        assert(var_size);
        double varValue = pow(0.5, rank * CONSTANT_FOR_RANK_CALC / var_size); //この式はSSIの定義次第で変更すること
        double polarity = phases[var_index];

        array<double, 3> element = {rank, polarity, varValue};
        csd[var_index] = element;
    }

    return csd;
}

double _get_non_zero_var_size_in_CSD(vector<array<double, 3>> csd)
{
    int counter = 0;
    for (int i = 0; i < (int)csd.size(); i++)
    {
        if (csd[i][0] != 0)
            counter++;
    }
    return (double)counter;
}

double calculate_SSI(vector<array<double, 3>> csd1, vector<array<double, 3>> csd2)
{
    if (csd1.size() != csd2.size() || csd1.size() == 0 || csd2.size() == 0)
        return 0;

    double ssi = 0;
    double size1 = _get_non_zero_var_size_in_CSD(csd1);
    double size2 = _get_non_zero_var_size_in_CSD(csd2);

    for (int i = 0; i < (int)csd1.size(); i++)
    {
        array<double, 3> val1 = csd1[i];
        array<double, 3> val2 = csd2[i];
        if (val1[0] == 0 || val2[0] == 0)
            continue;

        double rank1 = val1[0];
        double rank2 = val2[0];
        double phase1 = val1[1];
        double phase2 = val2[1];
        double value1 = val1[2];
        double value2 = val2[2];

        double similarity = (1 - abs(rank1 / size1 - rank2 / size2)) * (phase1 == phase2);
        double importance = 1 - abs(value1 - value2);
        ssi += similarity * importance;
    }

    if (size1 != 0 && size2 != 0)
    {
        double min_size = min(size1, size2);
        ssi /= min_size; //ノーマライゼーション
    }
    return ssi;
}

void _save_SSI(double ssi)
{
#pragma omp critical(SSI_DB)
    {
        SSI_database.push_back(ssi);
        if ((int)SSI_database.size() > LIMIT_SAVING_SSI)
            SSI_database.erase(SSI_database.begin());
    }
}
double _average(vector<double> v)
{
    double sum = 0;
    for (double s : v)
        sum += s;
    return sum / (double)v.size();
}
double _standardDeviation(vector<double> v)
{
    double sum2 = 0;
    for (double s : v)
        sum2 += s * s;
    double ave = _average(v);
    return sqrt(sum2 / (double)v.size() - ave * ave);
}

similarityLevel judge_SSI_score(double ssi)
{
    if (ssi == 0)
        return normal;

    vector<double> db;
#pragma omp critical(SHARED_CSD)
    db = SSI_database;

    double ave = _average(db);
    double std = _standardDeviation(db);
    _save_SSI(ssi); //次回以降の為に、今回のssiの値を保存

    if (ave == 0 || std == 0)
        return normal;
    if (ssi >= ave + std * ALPHA_TO_JUDGE_SSI || ((ave + std) >= 1 && ssi >= 0.99))
        return high;
    if (ssi < ave - std * ALPHA_TO_JUDGE_SSI || ((ave - std) <= 0 && ssi <= 0.01))
        return low;
    else
        return normal;
}

//UPDATE:: 同じ学習をした場合の2回目の実行時のファイル読み込み用
vector<int> restartNum;
vector<int> _split_string(const string &str, char delim = ',')
{
    istringstream iss(str);
    string tmp;
    vector<int> res;
    while (getline(iss, tmp, delim))
        res.push_back(stoi(tmp));
    return res;
}
vector<int> getSameLearntClause_restart(const string &path)
{
    stringstream ssurl{path};
    string buf;
    vector<string> strPath;
    while (std::getline(ssurl, buf, '/'))
        strPath.push_back(buf);

    string url = "/home/iida/experiments/20211102_sameLearntClause_similarityCheck/seed_2e5-targetRestart/" + strPath.back() + ".output";
    ifstream ifs(url);
    string line;
    while (getline(ifs, line))
        restartNum = _split_string(line);
    return restartNum;
}