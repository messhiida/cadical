#include "internal.hpp"

using namespace std;
vector<double> SSI_database;
vector<vector<array<double, 3>>> csd_database;
int sharedData = 0;
bool para_finished = false;

vector<vector<CaDiCaL::Clause *>> shared_learntClause(PARALLEL_NUM);
vector<vector<vector<array<double, 3>>>> shared_csd(PARALLEL_NUM);
//vector<bool> parallel_worker_action_table(PARALLEL_NUM, false);
array<array<bool, PARALLEL_NUM>, PARALLEL_NUM> parallel_worker_action_table;

bool check_action_table(int thread_num)
{
    int counter = 0;
    for (int i = 0; i < PARALLEL_NUM; i++)
    {
        if (parallel_worker_action_table[thread_num][i] == true)
            counter++;
    }
    if (counter >= ACTION_THREASHOLD)
        return true;
    else
        return false;
}

void set_bool_to_action_table(int thread_num, bool b)
{
    for (int i = 0; i < PARALLEL_NUM; i++)
    {
        parallel_worker_action_table[thread_num][i] = b;
        parallel_worker_action_table[i][thread_num] = b;
    }
}

void submit_csd(int thread_num, vector<array<double, 3>> csd)
{
    if ((int)csd.size() > 0)
        shared_csd[thread_num].push_back(csd);
    if ((int)shared_csd[thread_num].size() > LIMIT_SHARED_CSD)
        shared_csd[thread_num].erase(shared_csd[thread_num].begin());
}

//各workerごとにSSI table(自分, 自分以外)があり、それをもとにhigh / lowを判定
int update_worker_action_table()
{
    //clock_t start = clock();
    int counter = 0;
    for (int i = 0; i < PARALLEL_NUM; i++)
    {
        for (int j = 0; j < i; j++)
        {
            if (!shared_csd[i].empty() && !shared_csd[j].empty())
            {
                vector<array<double, 3>> csd1 = shared_csd[i].back();
                vector<array<double, 3>> csd2 = shared_csd[j].back();
                if (csd1.size() != csd2.size())
                    continue; //一旦サイズが異なればskipすることにしておく
                //if (csd1.size() == 0 || csd2.size() == 0)
                //printf("[%d][%d] %d, %d: %d, %d\n", i, j, csd1.size(), csd2.size(), shared_csd[i].back().size(), shared_csd[i].back().size());
                double ssi = calculate_SSI(csd1, csd2);
                similarityLevel sl = judge_SSI_score(ssi);
                if (sl == high && parallel_worker_action_table[i][j] == false)
                {
                    counter++;
                    printf("Similarity high: [%d][%d]\n", j, i);
                    parallel_worker_action_table[i][j] = true;
                    parallel_worker_action_table[j][i] = true;
                }
                if (sl == low && parallel_worker_action_table[i][j] == true)
                {
                    counter--;
                    printf("Similarity low: [%d][%d]\n", j, i);
                    parallel_worker_action_table[i][j] = false;
                    parallel_worker_action_table[j][i] = false;
                }
            }
        }
    }
    /*
    clock_t finish = clock();
    double spent = (double)(finish - start) / CLOCKS_PER_SEC;
    printf("Time spent: %lf, counter %d\n", spent, counter);
    */
    return counter;
}

vector<double> change_search_space(vector<double> score_table, double inc)
{
    double incrementalScore = CHANGE_SCORE_INCRE * inc; //score_incのINCRE回分をup (CIRと同じであれば10000回）
    while (incrementalScore > 1e150)                    //double上限には引っかからないとは思うが、念の為
        incrementalScore /= 10;

    for (double i = 0; i < (double)score_table.size() * CHANGE_SCORE_RATIO; i++)
    {
        //Minを探すことで、スコアを下から上に上げていっても上げていっても影響はないはず
        vector<double>::iterator minIt = min_element(score_table.begin(), score_table.end());
        size_t minIndex = distance(score_table.begin(), minIt);
        score_table[minIndex] += incrementalScore;
    }
    return score_table;
}

void submit_shared_learntClause(int thread_num, CaDiCaL::Clause *lc)
{
    for (int i = 0; i < PARALLEL_NUM; i++)
    {
        if (i != 0 && i != thread_num) //0 == master node, thread_num == 自分自身 の為この２つはskip
            shared_learntClause[i].push_back(lc);
    }
}

vector<CaDiCaL::Clause *> import_shared_learntClause()
{
    int my_thread = omp_get_thread_num();
    vector<CaDiCaL::Clause *> res;
    auto itr = shared_learntClause[my_thread].begin();
    while (itr != shared_learntClause[my_thread].end())
    {
        res.push_back(*itr);
        itr = shared_learntClause[my_thread].erase(itr);
    }
    return res;
}

int set_parallel_seed(int thread_num)
{
    int para_seed = 0;
    para_seed = 2048 * thread_num; //2048はmagic number
    return para_seed;
}

double _count_validScoreVars(vector<double> scores)
{
    double count = 0;
    for (double s : scores)
        if (s > CSD_SCORE_CRITERIA)
            count++;
    return count;
}

vector<array<double, 3>> get_CSD(vector<double> scores, vector<signed char> phases)
{
    int var_size = (int)scores.size();
    vector<array<double, 3>> csd(var_size); //csd[var] = {rank, phase, value},　不足分=下で定義されない分はzeroで埋められる
    vector<pair<double, int>> sortedScores; //ソートするための仮置場

    for (int i = 0; i < var_size; i++)
        sortedScores.push_back(make_pair(scores[i] * (double)(-1.0), i)); //ソートするためにスコアに(-1)をかける。後で反転して戻す
    sort(sortedScores.begin(), sortedScores.end());

    double size = _count_validScoreVars(scores);
    double rank = 0;

    for (int i = 0; i < (int)sortedScores.size(); i++)
    {
        pair<double, int> p = sortedScores[i];
        double score = p.first * (-1);
        int var_index = p.second;
        if (score <= CSD_SCORE_CRITERIA)
            break; //ソート済みの為、ここでbreakすることで一定以下のものはCSDの中に入らない
        rank++;
        double varValue = pow(0.5, rank * CONSTANT_FOR_RANK_CALC / size); //この式はSSIの定義次第で変更すること
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
    double ssi = 0;
    double size1 = _get_non_zero_var_size_in_CSD(csd1);
    double size2 = _get_non_zero_var_size_in_CSD(csd2);

    for (int i = 0; i < (int)csd1.size(); i++)
    {
        if (csd1.size() != csd2.size())
        {
            printf("size at i(var)=%d: csd1 %d, csd2 %d\n", i, (int)csd1.size(), (int)csd2.size());
            break; //skipする
        }
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

void save_CSD(vector<array<double, 3>> csd)
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
    if (ssi >= ave + std * ALPHA_TO_JUDGE_SSI || ((ave + std) >= 1 && ssi >= 0.99))
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