#ifndef PSF_MOMENTS_INCLUDED
#define PSF_MOMENTS_INCLUDED

#include <vif.hpp>
#include "filters.hpp"
#include "metrics.hpp"

using namespace vif;
using namespace vif::astro;

struct psf_options {
    vec1s libraries;
};

struct psf_library {
    vec1d lam, q11, q12, q22, w;
    filter_t vis_filter;
};

struct psf_moments {
    filter_database& filter_db;

    vec<1,psf_library> libraries;

    uint_t size() const {
        return libraries.size();
    }

    explicit psf_moments(filter_database& db) : filter_db(db) {}

    void read_options(program_arguments& opts) {
        vec1s psf_files;
        opts.read(arg_list(psf_files));

        psf_options popts;
        popts.libraries = psf_files;

        initialize(popts);
    }

    void initialize(const psf_options& opts) {
        filter_t vis_filter = filter_db.read_filter("euclid-vis");

        // Read monochromatic PSF libraries
        libraries.resize(opts.libraries.size());
        for (uint_t i : range(opts.libraries)) {
            auto& lib = libraries[i];

            fits::read_table(opts.libraries[i],
                "lambda", lib.lam, "w", lib.w, "q11", lib.q11, "q12", lib.q12, "q22", lib.q22
            );

            // Copy VIS filter
            lib.vis_filter = vis_filter;

            // Match it to the PSF vis_filter
            lib.w   = interpolate(lib.w,   lib.lam, lib.vis_filter.lam);
            lib.q11 = interpolate(lib.q11, lib.lam, lib.vis_filter.lam);
            lib.q12 = interpolate(lib.q12, lib.lam, lib.vis_filter.lam);
            lib.q22 = interpolate(lib.q22, lib.lam, lib.vis_filter.lam);
            lib.lam = lib.vis_filter.lam;

            // Ignore data outside of standard bandpass (450-950)
            {
                vec1u idi = where(lib.lam >= 0.450 && lib.lam <= 0.950);
                lib.w = lib.w[idi];
                lib.q11 = lib.q11[idi];
                lib.q12 = lib.q12[idi];
                lib.q22 = lib.q22[idi];
                lib.lam = lib.lam[idi];
                lib.vis_filter.res = lib.vis_filter.res[idi];
                lib.vis_filter.lam = lib.vis_filter.lam[idi];
            }

            // Include PSF weighting flux loss in PSF filter response
            lib.vis_filter.res *= lib.w;
            double norm = integrate(lib.vis_filter.lam, lib.vis_filter.res);
            vif_check(norm > 0 && is_finite(norm), "invalid PSF library ", opts.libraries[i]);
            lib.vis_filter.res /= norm;
        }
    }

    void get_moments(const vec1d& tlam, const vec1d& tsed,
        vec1d& q11, vec1d& q12, vec1d& q22, vec1d& ftot, vec1d& rlam) const {

        vif_check(!libraries.empty(), "uninitialized PSF moments");

        uint_t n = libraries.size();
        ftot.resize(n); q11.resize(n); q12.resize(n); q22.resize(n); rlam.resize(n);
        for (uint_t i : range(libraries)) {
            auto& lib = libraries[i];
            auto& f = lib.vis_filter;
            ftot[i] = sed2flux(f.lam, f.res,         tlam, tsed);
            q11[i]  = sed2flux(f.lam, f.res*lib.q11, tlam, tsed)/ftot[i];
            q12[i]  = sed2flux(f.lam, f.res*lib.q12, tlam, tsed)/ftot[i];
            q22[i]  = sed2flux(f.lam, f.res*lib.q22, tlam, tsed)/ftot[i];
            rlam[i] = sed2flux(f.lam, f.res*f.lam,   tlam, tsed)/ftot[i];

            if (ftot[i] == 0.0) {
                q11[i] = median(lib.q11);
                q12[i] = median(lib.q12);
                q22[i] = median(lib.q22);
                rlam[i] = median(f.lam);
            }

            vif_check(is_finite(ftot[i]), "invalid total flux for PSF ", i);
            vif_check(is_finite(q11[i]),  "invalid Q11 for PSF ", i);
            vif_check(is_finite(q12[i]),  "invalid Q12 for PSF ", i);
            vif_check(is_finite(q22[i]),  "invalid Q22 for PSF ", i);
        }
    }

    void get_moments_same_grid(const vec1d& tsed,
        vec1d& q11, vec1d& q12, vec1d& q22, vec1d& ftot, vec1d& rlam) const {

        vif_check(!libraries.empty(), "uninitialized PSF moments");

        uint_t n = libraries.size();
        ftot.resize(n); q11.resize(n); q12.resize(n); q22.resize(n); rlam.resize(n);
        for (uint_t i : range(libraries)) {
            auto& lib = libraries[i];
            auto& f = lib.vis_filter;

            vec1d s = tsed*f.res;

            ftot[i] = 0.0; q11[i] = 0.0; q12[i] = 0.0; q22[i] = 0.0; rlam[i] = 0.0;
            for (uint_t l : range(1, f.lam.size())) {
                ftot.safe[i] += (s.safe[l] + s.safe[l-1]);
                q11.safe[i]  += (s.safe[l]*lib.q11.safe[l] + s.safe[l-1]*lib.q11.safe[l-1]);
                q12.safe[i]  += (s.safe[l]*lib.q12.safe[l] + s.safe[l-1]*lib.q12.safe[l-1]);
                q22.safe[i]  += (s.safe[l]*lib.q22.safe[l] + s.safe[l-1]*lib.q22.safe[l-1]);
                rlam.safe[i] += (s.safe[l]*f.lam.safe[l]   + s.safe[l-1]*f.lam.safe[l-1]);
            }

            if (ftot[i] == 0.0) {
                q11[i] = median(lib.q11);
                q12[i] = median(lib.q12);
                q22[i] = median(lib.q22);
                rlam[i] = median(f.lam);
            } else {
                q11[i] /= ftot[i]; q12[i] /= ftot[i]; q22[i] /= ftot[i]; rlam[i] /= ftot[i];
            }

            vif_check(is_finite(ftot[i]), "invalid total flux for PSF ", i);
            vif_check(is_finite(q11[i]),  "invalid Q11 for PSF ", i);
            vif_check(is_finite(q12[i]),  "invalid Q12 for PSF ", i);
            vif_check(is_finite(q22[i]),  "invalid Q22 for PSF ", i);
        }
    }

    void get_moments(const vec1d& tlam, const vec1d& tsed,
        vec1d& q11, vec1d& q12, vec1d& q22) const {

        vec1d ftot, rlam;
        get_moments(tlam, tsed, q11, q12, q22, ftot, rlam);
    }

    void get_moments_same_grid(const vec1d& tsed, vec1d& q11, vec1d& q12, vec1d& q22) const {
        vec1d ftot, rlam;
        get_moments_same_grid(tsed, q11, q12, q22, ftot, rlam);
    }

    void get_moments(const vec1d& tlam, const vec1d& tsed, vec<1,metrics>& m) const {
        vec1d q11, q12, q22, ftot, rlam;
        get_moments(tlam, tsed, q11, q12, q22, ftot, rlam);

        m.resize(libraries.size());
        for (uint_t i : range(libraries)) {
            m[i] = metrics(q11[i], q12[i], q22[i], rlam[i], ftot[i]);
        };
    }

    void get_moments_same_grid(const vec1d& tsed, vec<1,metrics>& m) const {
        vec1d q11, q12, q22, ftot, rlam;
        get_moments_same_grid(tsed, q11, q12, q22, ftot, rlam);

        m.resize(libraries.size());
        for (uint_t i : range(libraries)) {
            m[i] = metrics(q11[i], q12[i], q22[i], rlam[i], ftot[i]);
        };
    }
};

#endif
