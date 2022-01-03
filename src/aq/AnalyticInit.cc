/*
 * (C) Copyright 2020-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "aq/AnalyticInit.h"

#include "aq/aq_geovals_interface.h"
#include "aq/GeoVals.h"
#include "aq/Locations.h"

namespace aq {

static oops::AnalyticInitMaker<ObsTraits, AnalyticInit> makerAnalytic1_("vortices");

// -----------------------------------------------------------------------------
AnalyticInit::AnalyticInit(const Parameters_ & options) : options_(options)
{ }
// -----------------------------------------------------------------------------
/*! \brief GeoVaLs Analytic Initialization
 *
 * \details This aq::AnalyticInit constructor was introduced for AQ in May, 2018 
 * for use with the interpolation test.
 *
 */
void AnalyticInit::fillGeoVaLs(const Locations & locs,
                               GeoVals & geovals) const {
  if (locs.comm().rank() == 0) {
    aq_geovals_analytic_init_f90(geovals.toFortran(), locs, options_.toConfiguration());
  }
}
// -----------------------------------------------------------------------------
}  // namespace aq
