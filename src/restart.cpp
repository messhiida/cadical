#include "internal.hpp"

namespace CaDiCaL
{

  // As observed by Chanseok Oh and implemented in MapleSAT solvers too,
  // various mostly satisfiable instances benefit from long quiet phases
  // with less or almost no restarts.  We implement this idea by prohibiting
  // the Glucose style restart scheme in a geometric fashion, which is very
  // similar to how originally restarts were scheduled in MiniSAT and earlier
  // solvers.  We start with say 1e3 = 1000 (opts.stabilizeinit) conflicts of
  // Glucose restarts.  Then in a "stabilizing" phase we disable these
  // until 1e4 = 2000 conflicts (if 'opts.stabilizefactor' is '200' percent)
  // have passed. After that we switch back to regular Glucose style restarts
  // until again 2 times more conflicts than the previous limit are reached.
  // Actually, in the latest version we still restarts during stabilization
  // but only in a reluctant doubling scheme with a rather high interval.

  bool Internal::stabilizing()
  {
    if (!opts.stabilize)
      return false;
    if (stable && opts.stabilizeonly)
      return true;
    if (stats.conflicts >= lim.stabilize)
    {
      report(stable ? ']' : '}');
      if (stable)
        STOP(stable);
      else
        STOP(unstable);
      stable = !stable;
      if (stable)
        stats.stabphases++;
      PHASE("stabilizing", stats.stabphases,
            "reached stabilization limit %" PRId64 " after %" PRId64 " conflicts",
            lim.stabilize, stats.conflicts);
      inc.stabilize *= opts.stabilizefactor * 1e-2;
      if (inc.stabilize > opts.stabilizemaxint)
        inc.stabilize = opts.stabilizemaxint;
      lim.stabilize = stats.conflicts + inc.stabilize;
      if (lim.stabilize <= stats.conflicts)
        lim.stabilize = stats.conflicts + 1;
      swap_averages();
      PHASE("stabilizing", stats.stabphases,
            "new stabilization limit %" PRId64 " at conflicts interval %" PRId64 "",
            lim.stabilize, inc.stabilize);
      report(stable ? '[' : '{');
      if (stable)
        START(stable);
      else
        START(unstable);
    }
    return stable;
  }

  // Restarts are scheduled by a variant of the Glucose scheme as presented in
  // our POS'15 paper using exponential moving averages.  There is a slow
  // moving average of the average recent glucose level of learned clauses as
  // well as a fast moving average of those glues.  If the end of a base
  // restart conflict interval has passed and the fast moving average is above
  // a certain margin over the slow moving average then we restart.

  bool Internal::restarting()
  {
    if (!opts.restart)
      return false;
    if ((size_t)level < assumptions.size() + 2)
      return false;
    if (stabilizing())
      return reluctant;
    if (stats.conflicts <= lim.restart)
      return false;
    double f = averages.current.glue.fast;
    double margin = (100.0 + opts.restartmargin) / 100.0;
    double s = averages.current.glue.slow, l = margin * s;
    LOG("EMA glue slow %.2f fast %.2f limit %.2f", s, f, l);
    return l <= f;
  }

  // This is Marijn's reuse trail idea.  Instead of always backtracking to the
  // top we figure out which decisions will be made again anyhow and only
  // backtrack to the level of the last such decision or to the top if no such
  // decision exists top (in which case we do not reuse any level).

  int Internal::reuse_trail()
  {
    if (!opts.restartreusetrail)
      return assumptions.size();
    int decision = next_decision_variable();
    assert(1 <= decision);
    int res = assumptions.size();
    if (use_scores())
    {
      while (res < level &&
             score_smaller(this)(decision, abs(control[res + 1].decision)))
        res++;
    }
    else
    {
      int64_t limit = bumped(decision);
      while (res < level && bumped(control[res + 1].decision) > limit)
        res++;
    }
    int reused = res - assumptions.size();
    if (reused > 0)
    {
      stats.reused++;
      stats.reusedlevels += reused;
      if (stable)
        stats.reusedstable++;
    }
    return res;
  }

  void Internal::restart()
  {
    START(restart);
    stats.restarts++;
    stats.restartlevels += level;
    if (stable)
      stats.restartstable++;
    LOG("restart %" PRId64 "", stats.restarts);
    backtrack(reuse_trail());

    lim.restart = stats.conflicts + opts.restartint;
    LOG("new restart limit at %" PRId64 " conflicts", lim.restart);

    //UPDATE:: CSD for every N restarts
    /*
  if(stats.restarts > 0 && stats.restarts % 100 == 0){
    map<int, vector<double> > new_csd = get_CSD (stab, phases.saved);
    if(csd_database.size() >= 2){
      map<int, vector<double> > old_csd = csd_database.back();
      double ssi = calculate_SSI(new_csd, old_csd);
      cout << stats.restarts << ",";
      similarityLevel sim = judge_SSI_score(ssi);
    }
    save_CSD(new_csd);
  }
  */

    //UPDATE:: restartごとにSSIを計算する
    /*
  if(csd_database.size() == LIMIT_SAVING_CSD){
    int index = LIMIT_SAVING_CSD - 2;
    double ssi = calculate_SSI(csd_database[0], csd_database[index]);
    std::cout << ssi << ", ["<< stats.restarts<< "]" << endl;
    map<int, vector<double> > latest_csd = csd_database[index+1];
    csd_database.clear();
    save_CSD(latest_csd);  
  }
  */

    //UPDATE:: clausesの内容を書き出し
    //その後、Restart時のCSDと同じ学習節があった場合の状態を書き出し
    //std::cout << "[Restart: " << stats.restarts << "], ";
    //std::cout << endl;

    //UPDATE:: CSD get function
    /*
  map<int, vector<double> > csd = get_CSD (stab, phases.saved);
  for (const auto& [key, value] : csd) {
    std::cout << "{" << key << ",";
    for(auto v: value) std::cout << v << ",";
    std::cout << "}, ";
  }
  std::cout << endl;
  */
    /*
  if(clauses.size() >= 10){
    for(int i=0; i<10; i++) read_learntClause(clauses[i]);
  }
  */
    //UPDATE:: import learnt clause for parallel
    /*
    if (PARALLEL_NUM > 1)
    {
      vector<Clause *> shared_clauses = import_shared_learntClause();
      //printf("[thread %d] total %d clauses + import %d clauses, learnt %d\n", omp_get_thread_num(), clauses.size(), shared_clauses.size(), stats.learned.clauses);
      for (const auto &c : shared_clauses)
      {

        
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
        

        clauses.push_back(c);
        stats.learned.literals += c->size;
        stats.learned.clauses++;
        //watch_clause(c);
        
      }
    }
    */
    if (PARALLEL_NUM > 1)
    {
      int my_thread = omp_get_thread_num();
      vector<array<double, 3>> my_csd = get_CSD(stab, phases.saved);
      if ((int)my_csd.size() > 0)
        submit_csd(my_thread, my_csd);

      bool change = check_action_table(my_thread);
      if (change == true)
      {
        stab = change_search_space(stab, score_inc);
        set_bool_to_action_table(my_thread, false);
      }
    }

    report('R', 2);
    STOP(restart);
  }
}
