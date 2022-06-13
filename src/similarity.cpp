#include "internal.hpp"

using namespace std;
vector<double> SSI_database;
vector<CSD> CSD_database;
int conflict_counter = 0;
CSD conflict_CSD;

vector<int> set_qtab(CaDiCaL::Queue queue, CaDiCaL::Links link)
{
    vector<int> queue_table;
    int next = queue.unassigned;
    while (next)
    {
        queue_table.push_back(next);
        next = link[next].prev;
    }
    return queue_table;
}

bool get_phase(int idx, bool target, CaDiCaL::Phases phases)
{
    // decide_phase関数(decide.cpp)を参照に変更、デフォルトのオプションでのみ動作

    const int initial_phase = 1;
    int phase = 0;
    if (!phase)
        phase = phases.forced[idx]; // TODO swap?
    if (!phase && target)
        phase = phases.target[idx];
    if (!phase)
        phase = phases.saved[idx];
    if (!phase)
        phase = initial_phase;

    return (bool)phase;
}

CSD get_CSD(vector<double> scoreTable, vector<int> queueTable, bool stable, CaDiCaL::Phases phases, int64_t conflicts)
{
    if (STABLE_ONLY_MODE)
        stable = true;

    int var_size = (int)scoreTable.size();
    CSD csd = init_csd(var_size);
    csd.conflicts = conflicts;

    if (stable) // score mode
    {
        vector<pair<double, int>> sortedScoreTable;
        for (int i = 0; i < var_size; i++)
        {
            double s = scoreTable[i];
            if (s > CSD_SET_CRITERIA)
                sortedScoreTable.push_back(make_pair(s * (-1), i)); //降順ソートするために、-1をかけておく
        }
        sort(sortedScoreTable.begin(), sortedScoreTable.end());
        int nonZeroVars = (int)sortedScoreTable.size();
        for (int j = 0; j < nonZeroVars; j++)
        {
            csd_element e;
            // double score = sortedScoreTable[j].first * (-1);
            int var = sortedScoreTable[j].second;
            e.rank = (int)j + 1; // rankの最上位は0始まりの為、+1で補正
            e.phase = get_phase(var, stable, phases);
            // e.phase = phases.saved[var];
            csd.data[var] = e;
        }
        csd.nonZeroVars = nonZeroVars;
    }
    else // queue mode
    {
        int size = (int)queueTable.size();
        for (int i = 0; i < size; i++)
        {
            csd_element e;
            int var = queueTable[i];
            e.rank = i + 1;
            e.phase = get_phase(var, stable, phases);
            csd.data[var] = e;
        }
        csd.nonZeroVars = size;
    }

    return csd;
}

CSD init_csd(size_t var_size = 0)
{
    CSD csd;
    csd.nonZeroVars = 0;
    csd.conflicts = 0;
    csd_element init_e;
    init_e.phase = false;
    init_e.rank = 0;
    init_e.value = 0;
    csd.data.resize((int)var_size);
    fill(csd.data.begin(), csd.data.end(), init_e);
    return csd;
}

void save_CSD(CSD csd)
{
    CSD_database.push_back(csd);
    if ((int)CSD_database.size() > LIMIT_SAVING_CSD + 1)
        CSD_database.erase(CSD_database.begin());
}

CSD get_prevCSD(int i)
{
    int size = (int)CSD_database.size();
    if (size - i - 1 < 0)
        return init_csd();
    else
        return CSD_database[(size - i - 1)];
}

double calculate_SSI(CSD csd1, CSD csd2)
{
    double size1 = (double)csd1.nonZeroVars;
    double size2 = (double)csd2.nonZeroVars;
    if (size1 == 0 || size2 == 0) //何もCSDの中に入っていない場合に相当、エラー処理
        return 0;

    double ssi = 0;
    double normalization = 0;
    for (size_t i = 0; i < csd1.data.size(); i++)
    {
        csd_element val1 = csd1.data[i];
        csd_element val2 = csd2.data[i];

        if (val1.rank == 0 || val2.rank == 0) //修論に合わせるため、どちらかのCSDがVSIDSを持っていなければ対象外とする
            continue;

        double similarity = 0.0;
        similarity = (1 - abs((double)val1.rank / size1 - (double)val2.rank / size2)) * (val1.phase == val2.phase);

        val1.value = pow(0.5, (double)val1.rank * CONSTANT_FOR_RANK_CALC / size1);
        val2.value = pow(0.5, (double)val2.rank * CONSTANT_FOR_RANK_CALC / size2);

        double importance = (val1.value + val2.value) / 2.0;
        // if (similarity != 1.0)
        //     printf("calc SSI - [%d] %lf,%lf by %d,%d : %d,%d\n", i, similarity, importance, val1.rank, size1, val2.rank, size2);
        ssi += similarity * importance;
        normalization += importance;
    }

    if (normalization == 0.0)
        return 0;

    ssi /= normalization;
    return ssi;
}

vector<double> change_search_space(vector<double> &score_table, CaDiCaL::ScoreSchedule &scores, double inc)
{
    double incrementalScore = CHANGE_SCORE_INCRE * inc; // score_incのINCRE回分をup (CIRと同じであれば10000回）
    while (incrementalScore > 1e150)                    // double上限には引っかからないとは思うが、念の為
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

void _save_SSI(double ssi)
{
    SSI_database.push_back(ssi);
    if ((int)SSI_database.size() > LIMIT_SAVING_SSI + 1)
        SSI_database.erase(SSI_database.begin());
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

SimilarityLevel judge_SSI_score(double ssi)
{
    if (ssi == 0)
        return normal;

    vector<double> db;
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