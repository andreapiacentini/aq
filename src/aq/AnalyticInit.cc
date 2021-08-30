/*
 * (C) Copyright 2020-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "model/AnalyticInit.h"

#include "model/GomAQ.h"
#include "model/LocationsAQ.h"

namespace aq {

// -----------------------------------------------------------------------------
AnalyticInit::AnalyticInit(const eckit::Configuration & config): config_(config)
{ }
// -----------------------------------------------------------------------------
/*! \brief GeoVaLs Analytic Initialization
 *
 * \details This aq::AnalyticInit constructor was introduced for QG in May, 2018 
 * for use with the interpolation test.
 *
 */
void AnalyticInit::fillGeoVaLs(const LocationsAQ & locs,
                               GomAQ & geovals) const {
  // Optionally replace values with analytic init
  if (config_.has("analytic_init"))
    aq_gom_analytic_init_f90(geovals.toFortran(), locs, config_);
}
// -----------------------------------------------------------------------------
}  // namespace aq
