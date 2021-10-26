#include <iostream>
#include <map>
#include "internal.hpp"

using namespace std;
vector<double> SSI_database;
vector< map<int, vector<double> > > csd_database;

/*** パラメーター for SSI ***/
//VSIDSの値がx未満であれば価値がないと判断する。0=VSIDSに値が入っていればOK。1=直近の学習説に含まれているものに限定
int scoreCriteria = 1; 
double constantForRankCalc = 10; //パラメーター定数   
int limitOfSavingSSI = 100; //SSIを最大いくつまで保存しAve/標準偏差の計算に用いるか
int limitOfSavingCSD = 10; //CSDを最大いくつまで保存しておくか into CSD_database
double alphaToJudgeSSI = 1; //SSIの類似度が高い場合に探索を変更するパラメーター：1=100%x標準偏差 (1回分の標準偏差のみ)
double affectRatioForChangingSearch = 0.1; //下位何％をpopするか
double incrementalNumForChangingSearch = 10000; //10000回分のVSIDS bumpに相当する分の変更を加える


void showResult(vector<double> s){
    //cout << "hello World :" << testCount /*<< stab.size()*/ << endl;
    for (int i=0; i < (int) s.size(); i++){
      cout << s[i] << ",";
    }
    cout << endl;
}

void draw_VSIDS_graph (vector<double> s){
    int varSize = (int) s.size();
    int onePlusVar = 0; //VSIDS scoreが1以上のvarを格納
    for (int i=0; i< (int) varSize; i++){
        if(s[i] >=1)    onePlusVar++;
    }
    cout << (double) onePlusVar / (double) varSize << ",";
}

double countValidScoreVariables(vector<double> scores){
    double count = 0;
    for(double s: scores) if(s > scoreCriteria) count++;
    return count;
}

map<int, vector<double> > get_CSD (vector<double> scores, vector<signed char> phases){
    map<int, vector<double> > csd; //csd[var] = {rank, phase, value}

    double size = countValidScoreVariables(scores);
    double rank = 0;
    
    //デバッグ用print
    //cout << "Total Size:" << scores.size() << ", valid:" << size << endl;    
    for(;;){ //無限ループだが、最大の値をscoreCriteriaに順々に書き換えていくことでbreakする
        vector<double>::iterator maxIt = max_element(scores.begin(), scores.end());
        
        if(*maxIt <= scoreCriteria) break; //全ての値をとった（入れ替えた）あとになっている為、ここでBreak
        
        rank++;
        double varValue = pow(0.5, rank * constantForRankCalc / size); //この式はSSIの定義次第で変更すること

        size_t maxIndex = distance(scores.begin(), maxIt);
        double polarity = phases[maxIndex];
        csd[maxIndex] = {(double) rank, polarity, varValue};
        
        //デバッグ用print
        //if (rank == 10 || rank == 100 || rank == 200) cout << "[" << rank << "]\t" << maxIndex << ":\t" << rank << ",\t" << (int) polarity << ",\t" << varValue << ",\t" << size << ",\t" << *maxIt << endl;
        
        *maxIt = scoreCriteria; //無限ループを避ける為SCORECRITERIALの値に入れ替え
    }
    
    return csd; 
}

double calculate_SSI (map<int, vector<double> > csd1, map<int, vector<double> > csd2){
    
    double SSI = 0;
    double size1 = csd1.size(); double size2 = csd2.size();

    for (const auto& [key1, val1]: csd1){
        auto itr = csd2.find(key1);
        if(itr != csd2.end()){
            vector<double> val2 = itr->second;
        
            double rank1 = val1[0]; double rank2 = val2[0]; 
            double phase1 = val1[1]; double phase2 = val2[1];
            double value1 = val1[2]; double value2 = val2[2];            

            double similarity = (1-abs(rank1/size1 - rank2/size2))*(phase1==phase2);
            double importance = abs(value1 - value2 - 1);
            SSI += similarity * importance;

            //デバッグ用print
            /*
            int key2 = itr->first;
            if(key1 == key2 && (key1 <= 10)){
                cout << key1 << ":" << key2 << ":\t";
                //cout << rank1 << "," << rank2 << ":\t" << phase1 << "," << phase2 << ":\t" << value1 << "," << value2 << ":\t" << size1 << "," << size2 << endl;
            } 
            */           
        }
    }
    
    if(csd1.size() != 0) SSI /= csd1.size(); //ノーマライゼーション: csd.size = 1x変数の数 = Maxの値のはずのため
    return SSI;
}

void save_SSI (double ssi){
    SSI_database.push_back(ssi);
    if((int) SSI_database.size() > limitOfSavingSSI) SSI_database.erase(SSI_database.begin());
}
void save_CSD (map<int, vector<double> > csd){
    if (csd.size() != 0) csd_database.push_back(csd); //最初の方はCSDが適切に取れず0になる為ifで例外処理
    if((int) csd_database.size() > limitOfSavingCSD) csd_database.erase(csd_database.begin());
}
double average(vector<double> v){
    double sum=0;
    for(double s: v) sum += s;
    return sum/(double)v.size();
}
double standardDeviation(vector<double> v){
    double sum2=0;
    for(double s: v) sum2 += s*s;
    double ave = average(v);
    return sqrt(sum2/(double)v.size() - ave*ave);
}
enum similarityLevel {high, normal, low};
similarityLevel judge_SSI_score(double ssi){

    double ave = average(SSI_database);
    double std = standardDeviation(SSI_database);
    save_SSI(ssi);

    similarityLevel res = normal; 
    if (ssi >= ave + std * alphaToJudgeSSI) res = high;
    if (ssi < ave - std * alphaToJudgeSSI) res = low;

    if ( (ave+std) >= 1 & ssi >= 0.99) res = high;
    if ( (ave-std) <= 0 && ssi <= 0.01) res = low;
    if ( ave == 0 || std == 0 ) res = normal;
    
    //デバッグ用print
    cout << "ave: " << ave << ", std: " << std << ", SSI: " << ssi << ", res: " << res << endl;
    return res;
}

//Todo ここまだこれから
/*
vector<double> change_search_space(vector<double> scores, double scoreInc){

    double incrementalScore = incrementalNumForChangingSearch * scoreInc; //score_inc 10000回分をup (CIRと同じ）

    for (double i=0; i < (double)scores.size()* (double) affectRatioForChangingSearch; i++){
        vector<double>::iterator minIt = min_element(scores.begin(), scores.end());
        size_t minIndex = distance(scores.begin(), minIt);
        scores[minIndex] += incrementalScore;
    }

    return scores;
}
*/
