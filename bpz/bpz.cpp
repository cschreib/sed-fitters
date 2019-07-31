#include "bpz.hpp"
#include "catalog_fitter.hpp"

int vif_main(int argc, char* argv[]) {
    filter_database db;
    psf_moments psf(db);
    bpz fitter(db, psf);
    catalog_fitter cat(fitter);

    // Read setup
    {
        program_arguments opts(argc, argv);

        // Override options
        {
            bool filter_flambda = true; // equivalent to FILTER_FORMAT=1
            bool filter_photons = true; // equivalent to FILTER_FORMAT=1
            bool trim_filters = true;
            opts.write(arg_list(filter_flambda, filter_photons, trim_filters));
        }

        // Default options
        {
            std::string filter_db = "/home/cschreib/data/fits/filter-db/db.dat";
            opts.read(arg_list(filter_db));
        }

        // Setup filter database
        db.read_options(opts);

        // Setup PSF moments
        psf.read_options(opts);

        // Setup input catalog
        cat.read_options(opts);

        // Setup fitter
        fitter.read_options(opts);
    }

    cat.fit();

    return 0;
}
