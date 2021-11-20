#include "internal.hpp"

using namespace std;
vector<double> SSI_database;
vector<map<int, vector<double>>> csd_database;

double _count_validScoreVars(vector<double> scores)
{
    double count = 0;
    for (double s : scores)
        if (s > CSD_SCORE_CRITERIA)
            count++;
    return count;
}

map<int, vector<double>> get_CSD(vector<double> scores, vector<signed char> phases)
{
    map<int, vector<double>> csd;           //csd[var] = {rank, phase, value}
    vector<pair<double, int>> sortedScores; //ソートするための仮置場

    for (int i = 0; i < (int)scores.size(); i++)
        sortedScores.push_back(make_pair(scores[i] * (-1), i)); //ソートするためにスコアに(-1)をかける。後で反転して戻す
    sort(sortedScores.begin(), sortedScores.end());

    double size = _count_validScoreVars(scores);
    double rank = 0;

    for (int i = 0; i < (int)sortedScores.size(); i++)
    {
        pair<double, int> p = sortedScores[i];
        double score = p.first * (-1);
        int varIndex = p.second;
        if (score <= CSD_SCORE_CRITERIA)
            break; //ソート済みの為、、ここでbreakすることで一定以下のものはCSDの中に入らない
        rank++;
        double varValue = pow(0.5, rank * CONSTANT_FOR_RANK_CALC / size); //この式はSSIの定義次第で変更すること
        double polarity = phases[varIndex];
        csd[varIndex] = {(double)rank, polarity, varValue};
    }

    return csd;
}

double calculate_SSI(map<int, vector<double>> csd1, map<int, vector<double>> csd2)
{
    double ssi = 0;
    double size1 = csd1.size();
    double size2 = csd2.size();

    for (const auto &[key1, val1] : csd1)
    {
        auto itr = csd2.find(key1);
        if (itr != csd2.end())
        {
            vector<double> val2 = itr->second;

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
    }

    if (csd1.size() != 0 && csd2.size() != 0)
    {
        double min_size = min(csd1.size(), csd2.size());
        ssi /= min_size; //ノーマライゼーション
    }
    return ssi;
}

void save_CSD(map<int, vector<double>> csd)
{
    if (csd.size() != 0)
        csd_database.push_back(csd); //CSDのサイズがゼロの場合は不要になるため例外処理
    if ((int)csd_database.size() > LIMIT_SAVING_CSD)
        csd_database.erase(csd_database.begin()); //古いCSDを頭から削除
}

void _save_SSI(double ssi)
{
    SSI_database.push_back(ssi);
    if ((int)SSI_database.size() > LIMIT_SAVING_SSI)
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

similarityLevel judge_SSI_score(double ssi)
{
    double ave = _average(SSI_database);
    double std = _standardDeviation(SSI_database);
    _save_SSI(ssi); //次回以降の為に、今回のssiの値を保存

    if (ave == 0 || std == 0)
        return normal;
    if (ssi >= ave + std * ALPHA_TO_JUDGE_SSI || ((ave + std) >= 1 & ssi >= 0.99))
        return high;
    if (ssi < ave - std * ALPHA_TO_JUDGE_SSI || ((ave - std) <= 0 && ssi <= 0.01))
        return low;
    else
        return normal;
}

vector<int> convert_learntClause_to_vector(CaDiCaL::Clause *c)
{
    vector<int> lc;
    for (int i = 0; i < c->size; i++)
        lc.push_back(c->literals[i]);
    return lc;
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

//Todo ここまだこれから
vector<double> change_search_space(vector<double> scores, double scoreInc)
{

    double affectRatioForChangingSearch = 0.1;                            //下位何％をpopするか
    double incrementalNumForChangingSearch = 10000;                       //10000回分のVSIDS bumpに相当する分の変更を加える
    double incrementalScore = incrementalNumForChangingSearch * scoreInc; //score_inc 10000回分をup (CIRと同じ）

    for (double i = 0; i < (double)scores.size() * (double)affectRatioForChangingSearch; i++)
    {
        vector<double>::iterator minIt = min_element(scores.begin(), scores.end());
        size_t minIndex = distance(scores.begin(), minIt);
        scores[minIndex] += incrementalScore;
    }

    return scores;
}