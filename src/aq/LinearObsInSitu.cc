/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "aq/LinearObsInSitu.h"

#include <vector>

#include "aq/aq_insitu_interface.h"
#include "aq/GeoVals.h"
#include "aq/ObsAuxControl.h"
#include "aq/ObsAuxIncrement.h"
#include "aq/ObsSpace.h"
#include "aq/ObsVec.h"
#include "eckit/config/Configuration.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

// -----------------------------------------------------------------------------
namespace aq {
// -----------------------------------------------------------------------------
static LinearObsOpMaker<LinearObsInSitu> makerInSituTL_("InSitu");
// -----------------------------------------------------------------------------

LinearObsInSitu::LinearObsInSitu(const ObsSpace &, const eckit::Configuration & config)
  : varin_(std::vector<std::string>{config.getString("mod var")})
{
  oops::Log::trace() << "LinearObsInSitu created" << std::endl;
}

// -----------------------------------------------------------------------------

void LinearObsInSitu::setTrajectory(const GeoVals &, const ObsAuxControl &) {}

// -----------------------------------------------------------------------------

void LinearObsInSitu::simulateObsTL(const GeoVals & geovals, ObsVec & ovec,
                                  const ObsAuxIncrement & bias) const {
  aq_insitu_equiv_tl_f90(geovals.toFortran(), ovec.toFortran(), bias.insitu());
}

// -----------------------------------------------------------------------------

void LinearObsInSitu::simulateObsAD(GeoVals & geovals, const ObsVec & ovec,
                                  ObsAuxIncrement & bias) const {
  aq_insitu_equiv_ad_f90(geovals.toFortran(), ovec.toFortran(), bias.insitu());
}

// -----------------------------------------------------------------------------

void LinearObsInSitu::print(std::ostream & os) const {
  os << "Linear in situ observation operator";
}

// -----------------------------------------------------------------------------

}  // namespace aq
