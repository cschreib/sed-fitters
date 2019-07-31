#ifndef COLOUR_CUBE_INCLUDED
#define COLOUR_CUBE_INCLUDED

#include "fitter_base.hpp"
#include "rebin.hpp"
#include "common.hpp"

struct fitter_options {
    std::string cube_file;
    double min_mag_err = 0.0;
    double lambda_step = 1e-3;
    std::string interp_method = "mcspline";
};

class colour_cube : public fitter_base {
public :
    // colour_cube PSF library
    vec<2,metrics> cube_psfs; // [npsf,nmodel]

    // Fit models
    uint_t nmodel = npos;
    vec1d cube_z;            // [nmodel]
    vec1d cube_prior;        // [nmodel]
    vec2d cube_flux;         // [nmodel,nlam]
    vec1d cube_lambda;       // [nlam]

    // Internal variables (used by all galaxies)
    vec2d tpl_flux;

    // Internal variables (used by each galaxy)
    struct workspace {
        vec1d pmodel;
        vec1d pz, pzc;

        vec1d wflux, weight, wmodel;
    };

    workspace global;

    // Config
    double rel_err = 0.0;
    double lambda_step = 1e-3; // [um]
    std::string interp_method = "mcspline";
    std::string cube_file;


    explicit colour_cube(filter_database& db, psf_moments& pm) : fitter_base(db, pm) {
        code_name = "colour_cube";
    }

    void do_read_options(program_arguments& opts) override {
        fitter_options fopts;

        // Behavior
        fopts.min_mag_err = 0.0;
        fopts.lambda_step = 1e-3;
        fopts.interp_method = "mcspline";

        // Library
        fopts.cube_file = "cube.fits";

        opts.read(arg_list(fopts.min_mag_err, fopts.cube_file, fopts.interp_method));

        initialize(fopts);
    }

    void initialize(const fitter_options& fopts) {
        nband = filters.size();

        // Options on the fit
        cube_file = fopts.cube_file;

        // Relative error on flux, sets minimum uncertainties
        rel_err = fopts.min_mag_err*(log(10.0)/2.5);

        // Wavelength step size for re-constructed SED (in microns)
        lambda_step = fopts.lambda_step;

        // Interpolation method to use when re-constructing high-resolution SED
        interp_method = fopts.interp_method;

        // Build colour_cube library
        make_sed_library();
    }

    void make_sed_library() {
        note("reading SED library");

        vec2d cube_bflux;
        vec1s cube_bands;

        fits::read_table(cube_file,
            "flux", cube_bflux, "bands", cube_bands,
            "fluxmb", cube_flux, "lambda_mb", cube_lambda, "z", cube_z, "weight", cube_prior
        );

        nmodel = cube_z.size();

        tpl_flux.resize(nmodel, nband);
        for (uint_t b : range(filters)) {
            uint_t ib = where_first(cube_bands == bands[b]);
            vif_check(ib != npos, "missing band '"+bands[b]+"' in cube");
            tpl_flux(_,b) = cube_bflux(_,ib);
        }
    }

    fit_result do_fit(const vec1d& ftot) override {
        fit_result fr;

        workspace local;
        if (multi_threaded) {
            reset_workspace(local);
        }

        workspace& w = (multi_threaded ? local : global);

        // Create weighted photometry
        for (uint_t l : range(nband)) {
            double bftot = ftot.safe[l];
            w.weight.safe[l] = 1.0/sqrt(phot_err2.safe[l] + sqr(bftot*rel_err));
            w.wflux.safe[l] = bftot*w.weight.safe[l];
        }

        // Find best SED and z, and store chi2 to build p(z,SED)
        uint_t iml = npos;
        fr.chi2 = finf;
        for (uint_t im : range(nmodel)) {
            // Create weighted model flux
            double wmm = 0.0;
            double wfm = 0.0;
            double wff = 0.0;
            for (uint_t l : range(nband)) {
                double wm = tpl_flux.safe(im,l)*w.weight.safe[l];
                double wf = w.wflux.safe[l];
                wfm += wf*wm;
                wff += sqr(wf);
                wmm += sqr(wm);
            }

            // Compute max likelihood scaling factor and chi2
            double scale = wfm/wmm;
            double tchi2 = wff - scale*wfm;

            // Compare and store
            w.pmodel.safe[im] = tchi2;

            if (tchi2 < fr.chi2) {
                fr.chi2 = tchi2;
                fr.z_obs = cube_z.safe[im];
                iml = im;
            }
        }

        // Create likelihood, apply prior, compute marginalized redshift
        fr.z_obsm = 0.0;
        double tprob = 0.0;
        for (uint_t im : range(nmodel)) {
            w.pmodel.safe[im] = exp(-0.5*(w.pmodel.safe[im] - fr.chi2))*cube_prior.safe[im];
            fr.z_obsm += w.pmodel.safe[im]*cube_z.safe[im];
            tprob += w.pmodel.safe[im];
        }

        w.pmodel /= tprob;
        fr.z_obsm /= tprob;

        // Store outputs

        // Compute PSFs
        if (!no_psf) {
            fr.psf_obs.resize(psfs.size());
            fr.psf_obsm.resize(psfs.size());
            for (uint_t ip : range(psfs)) {
                // Maximum likelihood
                fr.psf_obs[ip] = cube_psfs(ip,iml);

                // Marginalization
                metrics& tm = fr.psf_obsm[ip];
                for (uint_t im : range(nmodel)) {
                    tm += w.pmodel.safe[im]*cube_psfs(ip,im);
                }
            }
        }

        // Compute SEDs
        // TODO: save SED

        return fr;
    }

    void reset_workspace(workspace& w) {
        w.pmodel.resize(nmodel);

        w.wflux.resize(nband);
        w.weight.resize(nband);
        w.wmodel.resize(nband);
    }

    void prepare_fit() override {
        // Pre-compute colour_cube template fluxes and PSF moments
        bool compute_moments = false;
        if (!no_psf) {
            cube_psfs.resize(psfs.size(), nmodel);
            compute_moments = true;
        }

        if (compute_moments) {
            note("pre-computing model fluxes");
            vec1d olam = psfs.libraries[0].lam;

            auto pg = progress_start(nmodel);
            for (uint_t im : range(nmodel)) {
                // Interpolate to high res
                vec1d bsed = cube_flux(im,_);

                vif_check(count(!is_finite(bsed)) == 0, "SED ", im, " has invalid data");
                vif_check(count(bsed != 0.0) != 0, "SED ", im, " is zero");

                vec1d osed;
                if (interp_method == "cst") {
                    osed = rebin_cst(bsed, cube_lambda, olam);
                } else if (interp_method == "lin") {
                    osed = rebin_trapz(bsed, cube_lambda, olam);
                } else if (interp_method == "spline") {
                    osed = rebin_spline3(bsed, cube_lambda, olam);
                } else if (interp_method == "mcspline") {
                    osed = rebin_mcspline(bsed, cube_lambda, olam);
                }

                vif_check(count(!is_finite(osed)) == 0, "re-sampled SED ", im, " has invalid data");
                vif_check(count(osed != 0.0) != 0, "re-sampled SED ", im, " is zero");

                // Compute PSF moments
                vec<1,metrics> m;
                psfs.get_moments_same_grid(osed, m);

                vif_check(count(!is_finite(get_e1(m))) == 0, "PSF for SED ", im, " has invalid data");

                cube_psfs(_,im) = m;

                // TODO: compute SED

                progress(pg);
            }
        }

        // Create workspace
        if (!multi_threaded) {
            reset_workspace(global);
        }
    }
};

#endif
