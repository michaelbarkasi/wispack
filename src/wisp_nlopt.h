#ifndef WISPACK_WISP_NLOPT_H
#define WISPACK_WISP_NLOPT_H

/*
 * Minimal C++ convenience wrapper around NLopt's C API, routed entirely through nloptr's
 * R_GetCCallable-registered routines (nloptrAPI.h) instead of linking directly against a
 * system NLopt install. This is what lets wispack depend on the CRAN package `nloptr`
 * (LinkingTo) for NLopt instead of a raw `pkg-config`/hardcoded-path system library --
 * portable across the platforms CRAN tests, since nloptr ships its own NLopt build (via
 * CMake) when no suitable system one is found. See cran-comments.md.
 *
 * This is deliberately NOT a drop-in stand-in for upstream `nlopt.hpp`: it implements only
 * the slice of `nlopt::opt`'s interface this package actually calls (construct, set
 * min/max objective with a std::vector-based objective function, set_ftol_rel, set_maxeval,
 * set_stopval, optimize). nloptr's C API does not register `nlopt_set_munge` or
 * `nlopt_get_errmsg`, which upstream nlopt.hpp's own `opt` class needs unconditionally
 * (`set_munge` in its constructor, to support copying an `opt` object). Neither is needed
 * here: an `opt` is never copied (constructed, configured, optimized once, destroyed), and
 * on failure this throws a generic message instead of nlopt_get_errmsg's detailed one --
 * that message text is only ever used for an optional diagnostic print at the call sites in
 * wspc.cpp, never inspected programmatically.
 */

#include <nloptrAPI.h>
#include <vector>
#include <stdexcept>
#include <new>
#include <string>
#include <cstring>
#include <cmath>

namespace wisp_nlopt {

  enum algorithm { LD_LBFGS = NLOPT_LD_LBFGS };

  enum result {
    FAILURE          = NLOPT_FAILURE,
    INVALID_ARGS     = NLOPT_INVALID_ARGS,
    OUT_OF_MEMORY    = NLOPT_OUT_OF_MEMORY,
    ROUNDOFF_LIMITED = NLOPT_ROUNDOFF_LIMITED,
    FORCED_STOP      = NLOPT_FORCED_STOP,
    SUCCESS          = NLOPT_SUCCESS,
    STOPVAL_REACHED  = NLOPT_STOPVAL_REACHED,
    FTOL_REACHED     = NLOPT_FTOL_REACHED,
    XTOL_REACHED     = NLOPT_XTOL_REACHED,
    MAXEVAL_REACHED  = NLOPT_MAXEVAL_REACHED,
    MAXTIME_REACHED  = NLOPT_MAXTIME_REACHED
  };

  // Objective function signature, matching upstream nlopt::vfunc.
  typedef double (*vfunc)(const std::vector<double>& x, std::vector<double>& grad, void* data);

  class opt {
  public:
    opt(algorithm a, unsigned n)
      : o_(nlopt_create(static_cast<nlopt_algorithm>(a), n)), thunk_(NULL) {
      if (!o_) throw std::bad_alloc();
    }

    ~opt() {
      if (o_) nlopt_destroy(o_);
      delete thunk_;
    }

    opt(const opt&) = delete;
    opt& operator=(const opt&) = delete;

    void set_min_objective(vfunc f, void* data) { set_objective(f, data, true); }
    void set_max_objective(vfunc f, void* data) { set_objective(f, data, false); }

    void set_ftol_rel(double tol) { nlopt_set_ftol_rel(o_, tol); }
    void set_maxeval(int n)       { nlopt_set_maxeval(o_, n); }
    void set_stopval(double v)    { nlopt_set_stopval(o_, v); }

    // nloptr does not register nlopt_get_numevals as a C-callable, so this counts objective
    // evaluations itself (via the trampoline below) rather than asking NLopt for its own
    // count. Equivalent here since this package never sets any constraint (only an
    // objective plus bounds/tolerances), so the objective is the only thing NLopt evaluates.
    int get_numevals() const { return thunk_ ? thunk_->num_evals : 0; }

    result optimize(std::vector<double>& x, double& opt_f) {
      if (nlopt_get_dimension(o_) != x.size())
        throw std::invalid_argument("dimension mismatch");
      if (thunk_) thunk_->exception_pending = false;
      nlopt_result ret = nlopt_optimize(o_, x.empty() ? NULL : &x[0], &opt_f);
      if (thunk_ && thunk_->exception_pending) {
        thunk_->exception_pending = false;
        throw std::runtime_error(thunk_->exception_what);
      }
      if (ret < 0) {
        throw std::runtime_error(
          "NLopt optimization failed (code " + std::to_string(static_cast<int>(ret)) + ")"
        );
      }
      return static_cast<result>(ret);
    }

  private:
    struct thunk {
      vfunc f;
      void* data;
      nlopt_opt o;
      std::vector<double> xtmp, gradtmp, gradtmp0;
      bool exception_pending = false;
      std::string exception_what;
      int num_evals = 0;
    };

    nlopt_opt o_;
    thunk* thunk_;

    void set_objective(vfunc f, void* data, bool minimize) {
      delete thunk_;
      thunk_ = new thunk();
      thunk_->f = f;
      thunk_->data = data;
      thunk_->o = o_;
      unsigned n = nlopt_get_dimension(o_);
      thunk_->xtmp.assign(n, 0.0);
      thunk_->gradtmp.assign(n, 0.0);
      if (minimize) nlopt_set_min_objective(o_, &opt::trampoline, thunk_);
      else          nlopt_set_max_objective(o_, &opt::trampoline, thunk_);
    }

    // Raw nlopt_func callback: never let a C++ exception cross this boundary back into
    // NLopt's own (possibly non-C++) call stack. Catch it, force a graceful stop, and
    // re-throw from optimize() once nlopt_optimize() has actually returned.
    static double trampoline(unsigned n, const double* x, double* grad, void* data) {
      thunk* t = reinterpret_cast<thunk*>(data);
      t->num_evals++;
      try {
        if (n) std::memcpy(&t->xtmp[0], x, n * sizeof(double));
        double val = t->f(t->xtmp, grad ? t->gradtmp : t->gradtmp0, t->data);
        if (grad && n) std::memcpy(grad, &t->gradtmp[0], n * sizeof(double));
        return val;
      } catch (std::exception& e) {
        t->exception_pending = true;
        t->exception_what = e.what();
      } catch (...) {
        t->exception_pending = true;
        t->exception_what = "unknown exception in NLopt objective function";
      }
      nlopt_force_stop(t->o);
      return HUGE_VAL;
    }
  };

} // namespace wisp_nlopt

#endif
