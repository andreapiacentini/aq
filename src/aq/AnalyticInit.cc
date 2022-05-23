/*
 * This file is part of the Air Quality Ensemble Data Assimilation suite AQ.
 *
 * (C) Copyright 2022 CERFACS
 *
 * AQ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * any later version.
 *
 * AQ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * A copy of the GNU Lesser General Public License is distributed
 * along with AQ (files LICENSE.md, COPYING and COPYING.LESSER).
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
