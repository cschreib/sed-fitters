#ifndef FITTER_BASE_INCLUDED
#define FITTER_BASE_INCLUDED

#include <vif.hpp>
#include "filters.hpp"
#include "psf_moments.hpp"
#include "metrics.hpp"
#include "rebin.hpp"

using namespace vif;
using namespace vif::astro;

struct fit_result {
    float chi2;
    float z_obs, z_obsm;
    vec1f sed_obs, sed_obsm;
    vec<1,metrics> psf_obs, psf_obsm;
};

class fitter_base {
public :
    // Config
    filter_database& filter_db;
    psf_moments& psfs;
    std::string code_name;
    vec1s bands;
    vec1f phot_err2;
    uint_t nband = npos;
    vec<1,filter_t> filters;
    bool multi_threaded = false;
    bool no_psf = false;
    bool save_sed = false;
    std::string sed_interp_method = "cst";
    vec1d save_sed_lambda;

    explicit fitter_base(filter_database& db, psf_moments& pm) : filter_db(db), psfs(pm) {}

    virtual ~fitter_base() = default;

    void read_options(program_arguments& opts) {
        uint_t nthread = 0;
        vec1f depths;
        vec1s psf_files;
        double sed_lmin = 0.45;
        double sed_lmax = 0.95;
        double sed_step = 0.01;

        opts.read(arg_list(
            bands, nthread, psf_files, depths, save_sed, sed_lmin, sed_lmax, sed_step,
            sed_interp_method
        ));

        // Setup SED grid for saving SEDs
        if (save_sed) {
            uint_t npt = floor((sed_lmax-sed_lmin)/sed_step) + 1;
            save_sed_lambda = sed_step*indgen<double>(npt) + sed_lmin;

            vif_check(is_any_of(sed_interp_method, vec1s{"cst", "lin", "spline", "mcspline"}),
                "'sed_interp_method' must be one of 'cst', 'lin', 'spline', or 'mcspline' (got ",
                sed_interp_method, ")");
        }

        no_psf = psf_files.empty();

        multi_threaded = nthread != 0;

        // Find flux filters
        nband = bands.size();
        filters.resize(nband);
        for (uint_t l : range(bands)) {
            filters[l] = filter_db.read_filter(bands[l]);
        }

        // Square of photometric error (Gaussian additive component)
        phot_err2 = sqr(mag2uJy(depths)/10.0);
        vif_check(phot_err2.size() == nband, "mismatch between filters (", nband, ") and depths (",
            phot_err2.size(), ")");

        // Forward configuration to implementation
        do_read_options(opts);
    }

    virtual void do_read_options(program_arguments& opts) {}

    virtual void prepare_fit() = 0;

    virtual fit_result do_fit(const vec1d& ftot) = 0;

    fit_result do_fit(vec1d ftot, vec1d ferr) {
        // Check validity of input fluxes and correct accordingly
        vec1u idb = where(!is_finite(ftot) || !is_finite(ferr) || ferr < 0);
        ftot[idb] = 0;
        ferr[idb] = 1e9;

        // Set input error (overwrites the default)
        phot_err2 = sqr(ferr);

        return do_fit(ftot);
    }
};

#endif
