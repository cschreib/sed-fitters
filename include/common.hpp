#ifndef COMMON_INCLUDED
#define COMMON_INCLUDED

#include <vif.hpp>

vif::vec1d resample_sed(const vif::vec1d& tsed, const vif::vec1d& tlam, const vif::vec1d& blam) {
    double dl = blam[1] - blam[0];
    vif::uint_t npt = blam.size();
    vif::vec1d fobs(npt);
    for (uint_t l : vif::range(npt)) {
        if (blam.safe[l]-dl/2.0 > tlam.front() && blam.safe[l]+dl/2.0 < tlam.back()) {
            fobs.safe[l] = vif::integrate(tlam, tsed, blam.safe[l]-dl/2.0, blam.safe[l]+dl/2.0)/dl;
        }
    }

    return fobs;
}

#endif
