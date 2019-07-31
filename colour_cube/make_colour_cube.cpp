#include <vif.hpp>
#include "tree.hpp"

using namespace vif;
using namespace vif::astro;

int vif_main(int argc, char* argv[]) {
    vec1f ngal;
    vec1s bands_mb;
    vec1d lambda_mb;
    vec1f z, uv, vj;

    vec2f colors;
    vec2f fluxes;

    // Define the colors to use
    vec1s cbands = {"sdss-g","sdss-r","sdss-i","sdss-z","vista-Y","vista-J","vista-H"};

    // Define the width of color bins
    double dmag = 0.1;

    // Output file name
    std::string out = "color_cube.fits";
    std::string main_catalog;
    std::string mb_catalog;

    read_args(argc, argv, arg_list(cbands, dmag, main_catalog, mb_catalog, out));

    if (main_catalog.empty()) {
        error("must provide main_catalog=... with flux, redshift, and UVJ colors");
        return 1;
    }
    if (mb_catalog.empty()) {
        error("must provide mb_catalog=... with medium band fluxes");
        return 1;
    }

    vec1s cname;

    note("reading colors");

    vec1d norm;

    {
        vec1s bands;
        vec2f flux;
        fits::read_table(main_catalog, ftable(ngal, flux, bands, z, uv, vj));
        fits::read_table(mb_catalog, "bands", bands_mb, "lambda", lambda_mb);

        // Read in broad band fluxes
        norm.resize(flux.dims[0]);
        fluxes.resize(flux.dims[0], cbands.size());
        for (uint_t b : range(cbands)) {
            fluxes(_,b) = flux(_,where_first(bands == cbands[b]));
        }

        // Normalize to unit flux at reddest band
        for (uint_t i : range(ngal)) {
            norm[i] = fluxes(i,cbands.size()-1);
            fluxes(i,_) /= norm[i];
        }

        // Compute colors
        cname.resize(cbands.size()-1);
        colors.resize(flux.dims[0], cname.size());
        for (uint_t b : range(cbands.size()-1)) {
            cname[b] = cbands[b]+"/"+cbands[b+1];
            colors(_,b) = -2.5*log10(fluxes(_,b)/fluxes(_,b+1));
        }
    }

    uint_t ncol = cname.size();

    // Train
    note("design bins");

    vec<1,vec2f> bins(ncol);
    vec<1,vec1f> bcenter(ncol);
    vec<1,vec1f> bwidth(ncol);
    for (uint_t i : range(ncol)) {
        // Color is clamped to 3 to avoid long tails toward +inf for dropouts
        vec1f c = colors(_,i);
        c = clamp(c, -3.0, 3.0);

        // Select good data points
        vec1u ids = where(is_finite(c));

        // Get bounds
        double mi = min(c);
        double ma = max(c);

        // Define width of bins
        vec2f b = make_bins_from_edges(rgen_step(mi, ma, dmag));
        print(cname[i], ": ", mi, ", ", ma, ", ", b.dims[1]);

        bins[i] = b;
        bcenter[i] = bin_center(b);
        bwidth[i] = bin_width(b);
    }

    struct cell {
        double weight = 0.0;
        vec1d flux;
        vec1d fluxmb;
        double z = 0.0, zmin = +finf, zmax = -finf;
        double uv = 0.0, vj = 0.0;
    };

    // Setup averaging (using a sorted tree structure for fast lookup and low memory)
    tree<float,cell> mf;
    mf.setup(bins);

    fits::input_table itbl(mb_catalog);

    // Training data is loaded in chunks to avoid loading the entire data set in memory
    // at once, which could overflow the memory
    uint_t nchunk = 1e5;
    uint_t nstep = ceil(ngal.size()/float(nchunk));

    note("build tree");

    auto pg = progress_start(nstep);
    for (uint_t k : range(nstep)) {
        uint_t i0 = k*nchunk;
        uint_t i1 = (k+1)*nchunk;
        if (i1 > ngal.size()) i1 = ngal.size();
        --i1;

        vec1u cids = indgen(i1 - i0 + 1) + i0;

        // Grab fluxes from the file and store
        for (uint_t i : cids) {
            vec1f f;
            itbl.read_elements("flux", f, fits::at(i,_));
            f /= norm[i];

            // Remove invalid data (z=0 galaxies typically)
            if (!is_finite(f[0])) continue;

            // Find leaf of the tree with this color
            vec1u ik = mf.bin_key(colors(i,_));
            auto& c = mf.insert_binned(ik);

            // Average stuff using number density as weighting
            double w = ngal[i];
            c.weight += w;

            // Medium band fluxes
            if (c.fluxmb.empty()) {
                c.fluxmb = f*w;
            } else {
                c.fluxmb += f*w;
            }

            // Broad band fluxes
            if (c.flux.empty()) {
                c.flux = fluxes(i,_)*w;
            } else {
                c.flux += fluxes(i,_)*w;
            }

            // RF color
            c.uv += uv[i]*w;
            c.vj += vj[i]*w;

            // Redshift
            c.z += z[i]*w;
            c.zmin = min(c.zmin, z[i]);
            c.zmax = max(c.zmax, z[i]);
        }

        progress(pg);
    }

    note("normalize tree");

    // Finalize averaging, dividing by total weight
    uint_t nflx = 0;
    uint_t nflxmb = 0;
    mf.for_each_leaf([&](const vec1u& k, cell& c) {
        if (nflx == 0) nflx = c.flux.size();
        if (nflxmb == 0) nflxmb = c.fluxmb.size();

        c.flux /= c.weight;
        c.fluxmb /= c.weight;
        c.z /= c.weight;
        c.uv /= c.weight;
        c.vj /= c.weight;
    });

    uint_t ncell = mf.size();

    // Save binned data
    note("save tree (", ncell, " cells)");

    // First save meta data and allocate space
    {
        fits::output_table otbl(out);
        otbl.write_columns("colors", cname, "bands_mb", bands_mb, "lambda_mb", lambda_mb,
            "bands", cbands);

        for (uint_t i : range(ncol)) {
            otbl.write_column("bin"+to_string(i)+"_low", bins[i](0,_));
            otbl.write_column("bin"+to_string(i)+"_up", bins[i](1,_));
        }

        otbl.allocate_column<double>("flux",   ncell, nflx);
        otbl.allocate_column<double>("fluxmb", ncell, nflxmb);
        otbl.allocate_column<double>("weight", ncell);
        otbl.allocate_column<double>("z",      ncell);
        otbl.allocate_column<double>("uv",     ncell);
        otbl.allocate_column<double>("vj",     ncell);
        otbl.allocate_column<double>("zmin",   ncell);
        otbl.allocate_column<double>("zmax",   ncell);
        otbl.allocate_column<short>("id",      ncell, ncol);
    }

    // Now save individual cells in the sorted tree one by one
    fits::table otbl(out);
    uint_t il = 0;
    mf.for_each_leaf([&](const vec1u& k, cell& c) {
        otbl.update_elements("id",     k,        fits::at(il,_));
        otbl.update_elements("flux",   c.flux,   fits::at(il,_));
        otbl.update_elements("fluxmb", c.fluxmb, fits::at(il,_));
        otbl.update_elements("weight", c.weight, fits::at(il));
        otbl.update_elements("z",      c.z,      fits::at(il));
        otbl.update_elements("uv",     c.uv,     fits::at(il));
        otbl.update_elements("vj",     c.vj,     fits::at(il));
        otbl.update_elements("zmin",   c.zmin,   fits::at(il));
        otbl.update_elements("zmax",   c.zmax,   fits::at(il));
        ++il;
    });

    return 0;
}
