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

    if (PARALLEL_NUM > 1)
    {
      clock_t t1 = clock();

      int my_thread = omp_get_thread_num();
      vector<array<double, 3>> my_csd = get_CSD(stab, phases.saved, scores);

      clock_t t2 = clock();
      submit_csd(my_thread, my_csd);

      clock_t t3 = clock();
      bool change = check_ssi_table(my_thread);

      clock_t t4 = clock();
      if (change == true)
      {
        //TODO scoreがそもそもsortされていない件は解決していないが、大枠あっていそう（少なくともstabに入っているものは正しい）ので放置
        change_search_space(stab, scores, score_inc);
      }

      clock_t t5 = clock();

      printf("[%d]before: %d -> ", my_thread, (int)clauses.size());

      submit_shared_learntClause(my_thread, tmp_lc);
      vector<Clause *> tmp_slc = import_shared_learntClause(my_thread, clauses, stats);
      if (tmp_slc.size() != 0)
        for (const auto &c : tmp_slc)
        {
          if (proof)
            proof->add_derived_clause(c);
          watch_clause(c); //TODO ここがエラー
          if (c && opts.eagersubsume)
            eagerly_subsume_recently_learned_clauses(c);
        }

      printf("[%d]after: %d\n", my_thread, (int)clauses.size());

      if (stats.restarts % 100 == 0)
      {
        clock_t t6 = clock();
        double spent1 = (double)(t2 - t1) / CLOCKS_PER_SEC;
        double spent2 = (double)(t3 - t2) / CLOCKS_PER_SEC;
        double spent3 = (double)(t4 - t3) / CLOCKS_PER_SEC;
        double spent4 = (double)(t5 - t4) / CLOCKS_PER_SEC;
        double spent5 = (double)(t6 - t5) / CLOCKS_PER_SEC;
        printf("Time spent: [%d] %d, %lf, %lf, %lf, %lf, %lf, %d, %lld\n", my_thread, change, spent1, spent2, spent3, spent4, spent5, (int)stab.size(), stats.restarts);
      }
    }

    report('R', 2);
    STOP(restart);
  }
}
