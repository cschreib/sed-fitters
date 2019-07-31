#ifndef CATALOG_FITTER_INCLUDED
#define CATALOG_FITTER_INCLUDED

#include "fitter_base.hpp"
#include "common.hpp"

class catalog_fitter {
public :
    // Input SEDs
    uint_t ngal = 0;
    bool constant_error = false;

    // Config
    uint_t nthread = 0;
    fitter_base& fitter;
    uint_t npsf = 0;
    bool save_sed = false;
    uint_t sed_npt = 0;
    double min_mag_err = 0.01;
    std::string catalog;
    std::string out;
    vec1u filter_ids;

    explicit catalog_fitter(fitter_base& f) : fitter(f) {}

    void read_options(program_arguments& opts) {
        // Common fitter defaults (overrides the defaults in each fitter)
        min_mag_err = 0.01;

        // Read custom configuration from command line
        vec1s bands;
        vec1f depths;
        vec1s psf_files;
        opts.read(arg_list(
            min_mag_err, catalog, out, bands, psf_files, depths, nthread, save_sed, sed_npt
        ));

        npsf = psf_files.size();

        vif_check(!catalog.empty(), "please provide flux catalog in catalog=...");

        // Check input catalog
        note("checking input catalog data");
        fits::input_table itbl(catalog);

        // Check we have all the bands
        vec1s cat_bands;
        itbl.read_columns("bands", cat_bands);

        if (!bands.empty()) {
            filter_ids.resize(bands.size());
            for (uint_t i : range(bands)) {
                filter_ids[i] = where_first(cat_bands == bands[i]);
                vif_check(filter_ids[i] != npos, "could not find band '",
                    bands[i], "' in catalog");
            }

            cat_bands = bands;
        }

        // Get number of galaxies in catalog
        fits::column_info ci;
        vif_check(itbl.read_column_info("flux", ci), "missing 'FLUX' column");
        ngal = ci.dims[0];

        // Check if flux errors are provided
        bool had_no_depth = depths.empty();
        if (depths.empty() && itbl.read_column_info("depths", ci)) {
            itbl.read_column("depths", depths);
            depths = depths[filter_ids];
        }

        if (depths.empty()) {
            vif_check(itbl.read_column_info("flux_err", ci),
                "missing 'FLUX_ERR' or 'DEPTHS' columns (please set depths=... manually)");
        } else {
            vif_check(depths.size() == cat_bands.size(), "mismatch in depths and bands");
            fitter.phot_err2 = sqr(mag2uJy(depths)/10.0);
            constant_error = true;
        }

        if (out.empty()) {
            out = file::remove_extension(catalog)+"-zfit-"+fitter.code_name+".fits";
        }

        // Overwrite config options for fitter
        {
            if (had_no_depth) {
                if (depths.empty()) {
                    depths = replicate(99, cat_bands.size());
                }
                opts.write(arg_list(depths));
            }

            if (bands.empty()) {
                opts.write(arg_list(cat_bands));
            }
        }
    }

    void fit() {
        // Prepare fitter
        fitter.prepare_fit();

        // Open input catalog
        fits::input_table itbl(catalog);

        // Create output catalog
        note("creating output catalog");
        {
            file::mkdir(file::get_directory(out));
            fits::output_table otbl(out);

            vec1s id;
            if (itbl.read_column("id", id)) {
                otbl.write_column("id", id);
            }

            otbl.allocate_column<float>("chi2_obs", ngal);
            otbl.allocate_column<float>("z_obs",    ngal);
            otbl.allocate_column<float>("z_obsm",   ngal);

            if (save_sed) {
                otbl.allocate_column<float>("sed_obs",  ngal, sed_npt);
                otbl.allocate_column<float>("sed_obsm", ngal, sed_npt);
            }

            if (npsf != 0) {
                otbl.allocate_column<float>("e1_obs",    ngal, npsf);
                otbl.allocate_column<float>("e2_obs",    ngal, npsf);
                otbl.allocate_column<float>("r2_obs",    ngal, npsf);
                otbl.allocate_column<float>("rlam_obs",  ngal, npsf);
                otbl.allocate_column<float>("e1_obsm",   ngal, npsf);
                otbl.allocate_column<float>("e2_obsm",   ngal, npsf);
                otbl.allocate_column<float>("r2_obsm",   ngal, npsf);
                otbl.allocate_column<float>("rlam_obsm", ngal, npsf);
            }
        }

        fits::table otbl(out);
        std::mutex read_mutex;
        std::mutex write_mutex;

        // Function to fit one source
        auto do_source = [&](uint_t iter) {
            vec1d flux, flux_err;

            // Read data from input catalog
            {
                std::unique_lock<std::mutex> l(read_mutex);

                itbl.read_elements("flux", flux, fits::at(iter,_));
                if (!filter_ids.empty()) {
                    flux = flux[filter_ids];
                }

                if (!constant_error) {
                    itbl.read_elements("flux_err", flux_err, fits::at(iter,_));
                    if (!filter_ids.empty()) {
                        flux_err = flux_err[filter_ids];
                    }
                }
            }

            // Fit
            fit_result fr = (constant_error ? fitter.do_fit(flux) : fitter.do_fit(flux, flux_err));

            // Update output catalog
            {
                std::unique_lock<std::mutex> l(write_mutex);

                otbl.update_elements("chi2_obs", fr.chi2, fits::at(iter));
                otbl.update_elements("z_obs", fr.z_obs,   fits::at(iter));
                otbl.update_elements("z_obsm", fr.z_obsm, fits::at(iter));

                if (save_sed) {
                    otbl.update_elements("sed_obs",  fr.sed_obs,  fits::at(iter,_));
                    otbl.update_elements("sed_obsm", fr.sed_obsm, fits::at(iter,_));
                }

                if (npsf != 0) {
                    otbl.update_elements("e1_obs",    get_e1(fr.psf_obs),    fits::at(iter,_));
                    otbl.update_elements("e2_obs",    get_e2(fr.psf_obs),    fits::at(iter,_));
                    otbl.update_elements("r2_obs",    get_r2(fr.psf_obs),    fits::at(iter,_));
                    otbl.update_elements("rlam_obs",  get_rlam(fr.psf_obs),  fits::at(iter,_));
                    otbl.update_elements("e1_obsm",   get_e1(fr.psf_obsm),   fits::at(iter,_));
                    otbl.update_elements("e2_obsm",   get_e2(fr.psf_obsm),   fits::at(iter,_));
                    otbl.update_elements("r2_obsm",   get_r2(fr.psf_obsm),   fits::at(iter,_));
                    otbl.update_elements("rlam_obsm", get_rlam(fr.psf_obsm), fits::at(iter,_));
                }
            }
        };

        note("started fitting");
        thread::parallel_for pfor(nthread);
        pfor.verbose = true;
        pfor.progress_step = 1;
        pfor.update_rate = 2.0;
        pfor.chunk_size = 100;
        pfor.execute(do_source, ngal);
        note("fitting done");
    }
};

#endif
