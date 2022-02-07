#include "internal.hpp"

using namespace std;
vector<double> SSI_database;

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