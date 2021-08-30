/*
 * (C) Copyright 2019-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include<memory>

#include "aq/GetValuesTLAD.h"

#include "oops/util/Logger.h"

#include "aq/GeometryAQ.h"
#include "aq/GomAQ.h"
#include "aq/IncrementAQ.h"
#include "aq/LocationsAQ.h"
#include "aq/StateAQ.h"

namespace aq {

// -----------------------------------------------------------------------------
/// Constructor, destructor
// -----------------------------------------------------------------------------
GetValuesTLAD::GetValuesTLAD(const GeometryAQ & geom, const LocationsAQ & locs,
                             const eckit::Configuration & linearGetValuesConf)
  : locs_(locs)
{
  aq_getvalues_setup_f90(hmat_);
  oops::Log::trace() << "GetValuesTLAD create: linearGetValuesConf = " <<
                        linearGetValuesConf << std::endl;
}
// -----------------------------------------------------------------------------
GetValuesTLAD::~GetValuesTLAD() {
  aq_getvalues_delete_f90(hmat_);
}
// -----------------------------------------------------------------------------
void GetValuesTLAD::setTrajectory(const StateAQ & state, const util::DateTime & t1,
                                  const util::DateTime & t2, GomAQ & geovals) {
  aq_getvalues_build_f90(locs_, state.fields().toFortran(),
                         t1, t2, hmat_);
}
// -----------------------------------------------------------------------------
/// Get increment values at observation locations
// -----------------------------------------------------------------------------
void GetValuesTLAD::fillGeoVaLsTL(const IncrementAQ & inc, const util::DateTime & t1,
                                  const util::DateTime & t2, GomAQ & geovals) const {
  aq_getvalues_interp_tl_f90(locs_, inc.fields().toFortran(),
                             t1, t2, hmat_, geovals.toFortran());
}
// -----------------------------------------------------------------------------
void GetValuesTLAD::fillGeoVaLsAD(IncrementAQ & inc, const util::DateTime & t1,
                                  const util::DateTime & t2, const GomAQ & geovals) const {
  aq_getvalues_interp_ad_f90(locs_, inc.fields().toFortran(),
                             t1, t2, hmat_, geovals.toFortran());
}
// -----------------------------------------------------------------------------
void GetValuesTLAD::print(std::ostream & os) const {
  os << "AQ GetValues TL/AD";
}
// -----------------------------------------------------------------------------

}  // namespace aq
