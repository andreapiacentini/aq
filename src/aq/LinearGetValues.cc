/*
 * (C) Copyright 2019-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include<memory>

#include "aq/LinearGetValues.h"

#include "oops/util/Logger.h"

#include "aq/aq_getvalues_interface.h"
#include "aq/Geometry.h"
#include "aq/GeoVals.h"
#include "aq/Increment.h"
#include "aq/Locations.h"
#include "aq/State.h"

namespace aq {

// -----------------------------------------------------------------------------
/// Constructor, destructor
// -----------------------------------------------------------------------------
LinearGetValues::LinearGetValues(const Geometry & geom, const Locations & locs,
                             const eckit::Configuration & linearGetValuesConf)
  : locs_(locs)
{
  aq_getvalues_setup_f90(hmat_);
  oops::Log::trace() << "LinearGetValues create: linearGetValuesConf = " <<
                        linearGetValuesConf << std::endl;
}
// -----------------------------------------------------------------------------
LinearGetValues::~LinearGetValues() {
  aq_getvalues_delete_f90(hmat_);
}
// -----------------------------------------------------------------------------
void LinearGetValues::setTrajectory(const State & state, const util::DateTime & t1,
                                  const util::DateTime & t2, GeoVals & geovals) {
  aq_getvalues_build_f90(locs_, state.fields().toFortran(),
                         t1, t2, hmat_);
}
// -----------------------------------------------------------------------------
/// Get increment values at observation locations
// -----------------------------------------------------------------------------
void LinearGetValues::fillGeoVaLsTL(const Increment & inc, const util::DateTime & t1,
                                  const util::DateTime & t2, GeoVals & geovals) const {
  aq_getvalues_interp_tl_f90(locs_, inc.fields().toFortran(),
                             t1, t2, hmat_, geovals.toFortran());
}
// -----------------------------------------------------------------------------
void LinearGetValues::fillGeoVaLsAD(Increment & inc, const util::DateTime & t1,
                                  const util::DateTime & t2, const GeoVals & geovals) const {
  aq_getvalues_interp_ad_f90(locs_, inc.fields().toFortran(),
                             t1, t2, hmat_, geovals.toFortran());
}
// -----------------------------------------------------------------------------
void LinearGetValues::print(std::ostream & os) const {
  os << "AQ GetValues TL/AD";
}
// -----------------------------------------------------------------------------

}  // namespace aq
