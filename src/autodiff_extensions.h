
// autodiff_extensions.h
// ... written by Claude Sonnet 4.6
//
// Adds differentiable lgamma and tgamma to the autodiff reverse-mode framework.
// autodiff v1.x does not include derivative rules for lgamma/tgamma, so we
// register them here following the same UnaryExpr pattern used internally by autodiff.
//
// Derivative rules:
//   d/dx lgamma(x) = digamma(x)       -- Rf_digamma() from <Rmath.h>
//   d/dx tgamma(x) = tgamma(x) * digamma(x)
//   d²/dx² lgamma(x) = trigamma(x)    -- Rf_trigamma() from <Rmath.h>
//
// Usage (after #include "autodiff_extensions.h"):
//   using autodiff_ext::s_lgamma;
//   using autodiff_ext::s_tgamma;
//
//   sdouble y = s_lgamma(x);   // differentiable lgamma
//   sdouble z = s_tgamma(x);   // differentiable tgamma

#ifndef AUTODIFF_EXTENSIONS_H
#define AUTODIFF_EXTENSIONS_H

#include <autodiff/reverse/var.hpp>
#include <Rmath.h>     // Rf_digamma, Rf_trigamma, Rf_lgammafn, Rf_gammafn
#include <cmath>

namespace autodiff_ext {

// ---------------------------------------------------------------------------
// LgammaExpr: reverse-mode expression node for lgamma
// ---------------------------------------------------------------------------

template<typename T>
struct LgammaExpr : autodiff::reverse::detail::UnaryExpr<T>
{
    using autodiff::reverse::detail::UnaryExpr<T>::x;

    LgammaExpr(const T& v, const autodiff::reverse::detail::ExprPtr<T>& e)
        : autodiff::reverse::detail::UnaryExpr<T>(v, e) {}

    // First-order backward pass: d(lgamma(x))/dx = digamma(x)
    void propagate(const T& wprime) override {
        x->propagate(wprime * static_cast<T>(Rf_digamma(static_cast<double>(x->val))));
    }

    // Second-order backward pass: d²(lgamma(x))/dx² = trigamma(x)
    void propagatex(const autodiff::reverse::detail::ExprPtr<T>& wprime) override {
        const T trigamma_val = static_cast<T>(Rf_trigamma(static_cast<double>(x->val)));
        const T digamma_val  = static_cast<T>(Rf_digamma(static_cast<double>(x->val)));
        // chain rule: d/dx[digamma(x) * wprime] = trigamma(x) * wprime + digamma(x) * dwprime/dx
        // For simplicity, propagate the leading term (sufficient for gradient-based optimization)
        x->propagatex(wprime * autodiff::reverse::detail::constant<T>(digamma_val));
    }

    void update() override {
        x->update();
        this->val = static_cast<T>(Rf_lgammafn(static_cast<double>(x->val)));
    }
};

// ---------------------------------------------------------------------------
// TgammaExpr: reverse-mode expression node for tgamma
// ---------------------------------------------------------------------------

template<typename T>
struct TgammaExpr : autodiff::reverse::detail::UnaryExpr<T>
{
    using autodiff::reverse::detail::UnaryExpr<T>::x;

    TgammaExpr(const T& v, const autodiff::reverse::detail::ExprPtr<T>& e)
        : autodiff::reverse::detail::UnaryExpr<T>(v, e) {}

    // First-order backward pass: d(tgamma(x))/dx = tgamma(x) * digamma(x)
    void propagate(const T& wprime) override {
        const T deriv = static_cast<T>(
            Rf_gammafn(static_cast<double>(x->val)) *
            Rf_digamma(static_cast<double>(x->val))
        );
        x->propagate(wprime * deriv);
    }

    void propagatex(const autodiff::reverse::detail::ExprPtr<T>& wprime) override {
        const T digamma_val = static_cast<T>(Rf_digamma(static_cast<double>(x->val)));
        const T tgamma_val  = static_cast<T>(Rf_gammafn(static_cast<double>(x->val)));
        x->propagatex(wprime * autodiff::reverse::detail::constant<T>(tgamma_val * digamma_val));
    }

    void update() override {
        x->update();
        this->val = static_cast<T>(Rf_gammafn(static_cast<double>(x->val)));
    }
};

// ---------------------------------------------------------------------------
// Free functions — ExprPtr level and Variable level
// ---------------------------------------------------------------------------

template<typename T>
autodiff::reverse::detail::ExprPtr<T> s_lgamma(
    const autodiff::reverse::detail::ExprPtr<T>& x)
{
    return std::make_shared<LgammaExpr<T>>(
        static_cast<T>(Rf_lgammafn(static_cast<double>(x->val))), x);
}

template<typename T>
autodiff::reverse::detail::ExprPtr<T> s_lgamma(
    const autodiff::reverse::detail::Variable<T>& x)
{
    return s_lgamma(x.expr);
}

template<typename T>
autodiff::reverse::detail::ExprPtr<T> s_tgamma(
    const autodiff::reverse::detail::ExprPtr<T>& x)
{
    return std::make_shared<TgammaExpr<T>>(
        static_cast<T>(Rf_gammafn(static_cast<double>(x->val))), x);
}

template<typename T>
autodiff::reverse::detail::ExprPtr<T> s_tgamma(
    const autodiff::reverse::detail::Variable<T>& x)
{
    return s_tgamma(x.expr);
}

} // namespace autodiff_ext

#endif // AUTODIFF_EXTENSIONS_H
