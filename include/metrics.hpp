#ifndef METRICS_INCLUDED
#define METRICS_INCLUDED

#include <vif.hpp>

using namespace vif;

// Average PSF metrics
struct metrics {
    double q11 = 0, q12 = 0, q22 = 0;
    double e1 = 0, e2 = 0, r2 = 0;
    double rlam = 0, ftot = 0;

    metrics() = default;

    metrics(double) {} // to enable integrate()

    explicit metrics(double t11, double t12, double t22, double rl, double ft = 1.0) :
        q11(t11), q12(t12), q22(t22), rlam(rl), ftot(ft) {
        get_ellipticities();
    }

    void get_ellipticities() {
        r2 = q11+q22;
        e1 = (q11-q22)/r2;
        e2 = 2*q12/r2;
    }

    metrics& operator *= (double norm) {
        q11 *= norm; q12 *= norm; q22 *= norm;
        e1 *= norm; e2 *= norm; r2 *= norm;
        rlam *= norm; ftot *= norm;
        return *this;
    }

    metrics& operator *= (const metrics& m) {
        q11 *= m.q11; q12 *= m.q12; q22 *= m.q22;
        e1 *= m.e1; e2 *= m.e2; r2 *= m.r2;
        rlam *= m.rlam; ftot *= m.ftot;
        return *this;
    }

    metrics& operator /= (double norm) {
        q11 /= norm; q12 /= norm; q22 /= norm;
        e1 /= norm; e2 /= norm; r2 /= norm;
        rlam /= norm; ftot /= norm;
        return *this;
    }

    metrics& operator += (const metrics& m) {
        q11 += m.q11; q12 += m.q12; q22 += m.q22;
        e1 += m.e1; e2 += m.e2; r2 += m.r2;
        rlam += m.rlam; ftot += m.ftot;
        return *this;
    }

    metrics& operator -= (const metrics& m) {
        q11 -= m.q11; q12 -= m.q12; q22 -= m.q22;
        e1 -= m.e1; e2 -= m.e2; r2 -= m.r2;
        rlam -= m.rlam; ftot -= m.ftot;
        return *this;
    }

    void reset() {
        operator*=(0.0);
    }
};

metrics operator* (metrics m1, const metrics& m2) {
    return m1 *= m2;
}
metrics operator* (metrics m, double norm) {
    return m *= norm;
}
metrics operator* (double norm, metrics m) {
    return m *= norm;
}
metrics operator+ (metrics m1, const metrics& m2) {
    return m1 += m2;
}
metrics operator- (metrics m1, const metrics& m2) {
    return m1 -= m2;
}

metrics sqrt(metrics m) {
    m.q11 = sqrt(m.q11);
    m.q12 = sqrt(m.q12);
    m.q22 = sqrt(m.q22);
    m.e1 = sqrt(m.e1);
    m.e2 = sqrt(m.e2);
    m.r2 = sqrt(m.r2);
    m.rlam = sqrt(m.rlam);
    m.ftot = sqrt(m.ftot);
    return m;
}

namespace std {
    template<>
    struct is_arithmetic<metrics> : std::true_type {}; // to enable integrate()
}

auto get_q11 = vectorize_lambda([](const metrics& n) {
    return n.q11;
});
auto get_q12 = vectorize_lambda([](const metrics& n) {
    return n.q12;
});
auto get_q22 = vectorize_lambda([](const metrics& n) {
    return n.q22;
});
auto get_e1 = vectorize_lambda([](const metrics& n) {
    return n.e1;
});
auto get_e2 = vectorize_lambda([](const metrics& n) {
    return n.e2;
});
auto get_r2 = vectorize_lambda([](const metrics& n) {
    return n.r2;
});
auto get_rlam = vectorize_lambda([](const metrics& n) {
    return n.rlam;
});
auto get_ftot = vectorize_lambda([](const metrics& n) {
    return n.ftot;
});

#endif
