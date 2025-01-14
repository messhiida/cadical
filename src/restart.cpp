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

  // UPDATE:: 追加
  int luby(int x)
  {
    int size, seq;
    for (size = 1, seq = 0; size < x + 1; seq++, size = 2 * size + 1)
      ;
    while (size - 1 != x)
    {
      size = (size - 1) >> 1;
      seq--;
      x = x % size;
    }
    return pow(2, seq);
  }
  int lubyFunc(int restart)
  {
    int counter = 0;
    for (int i = 0; i < restart; i++)
      counter += luby(i);
    return counter;
  }
  int geometricFunc(int restart)
  {
    int counter = 0;
    for (int i = 0, repeat = 0; i < restart; i++, repeat++)
    {
      int c = (int)(100 * pow(1.5, repeat));
      if (c >= 3000)
      {
        repeat = 0;
        c = (int)(100 * pow(1.5, repeat));
      }
      counter += c;
      // printf("[counter] %d %d %d %d %d\n", i, repeat, c, counter, restart);
    }
    return counter;
  }

  bool Internal::restarting()
  {
    if (RESTART_POLICY == UNIFORM_INTERVAL)
    {
      int criteria = 256 * stats.restarts;
      if (stats.conflicts >= criteria)
        return true;
      else
        return false;
    }
    else if (RESTART_POLICY == GEOMETRIC_INTERVAL)
    {
      int criteria = geometricFunc(stats.restarts);
      // int criteria = 100 * pow(1.5, (stats.restarts - 1));
      if (stats.conflicts >= criteria)
      {
        // printf("%ld %d %ld\n", stats.conflicts, criteria, stats.restarts);
        return true;
      }
      else
        return false;
    }
    else if (RESTART_POLICY == LUBY_INTERVAL)
    {
      int criteria = 100 * lubyFunc(stats.restarts);
      if (stats.conflicts >= criteria)
      {
        // printf("%d conflict : %d %d %d\n", RESTART_POLICY, stats.conflicts, stats.restarts, criteria);
        return true;
      }
      else
        return false;
    }
    else
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

    // UPDATE::
    //ここからrestart時のSSI計算用
    /*
    vector<int> qtab = set_qtab(queue, links);
    CSD csd = get_CSD(stab, qtab, stable, phases, stats.conflicts);
    save_CSD(csd);
    if (stats.restarts % 10 == 0)
    {
      printf("SSI-r:[%d] %d ", (int)stats.restarts, stable);
      for (int j = 0; j < 3; j++)
      {
        for (int i = 1; i <= 10; i++)
        {
          int p = i * (int)pow(10, j);
          if (p >= (int)stats.restarts)
          {
            printf("- - ");
            continue;
          }
          CSD prev = get_prevCSD(p);
          double ssi = 0.0;
          ssi = calculate_SSI(csd, prev);
          int64_t confllicts_between = csd.conflicts - prev.conflicts;
          printf("%lf %lld ", ssi, confllicts_between);
        }
      }
      printf("\n");
    }
    */

    //ここからconflict時のSSI計算用
    // conflict_counter = 0;

    report('R', 2);
    STOP(restart);
  }
}
