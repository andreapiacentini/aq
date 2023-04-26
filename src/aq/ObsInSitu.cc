/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2021-2022 CERFACS.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "aq/ObsInSitu.h"

#include <vector>

#include "aq/aq_insitu_interface.h"
#include "aq/GeoVals.h"
#include "aq/Locations.h"
#include "aq/ObsAuxControl.h"
#include "aq/ObsSpace.h"
#include "aq/ObsVec.h"
#include "aq/TraitsFwd.h"

#include "eckit/config/Configuration.h"

#include "oops/base/Locations.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

// -----------------------------------------------------------------------------
namespace aq {
// -----------------------------------------------------------------------------
static ObsOpMaker<InSitu> makerInSitu_("InSitu");
// -----------------------------------------------------------------------------

InSitu::InSitu(const ObsSpace & odb, const eckit::Configuration & config)
  : obsdb_(odb), varin_(std::vector<std::string>{config.getString("mod var")})
{
  oops::Log::trace() << "InSitu created." << std::endl;
}

// -----------------------------------------------------------------------------

void InSitu::simulateObs(const GeoVals & geovals, ObsVec & ovec,
                              const ObsAuxControl & bias) const {
  if (obsdb_.comm().rank() == 0) {
    aq_insitu_equiv_f90(geovals.toFortran(), ovec.toFortran(), bias.insitu());
  }
}

// -----------------------------------------------------------------------------

ObsOpBase::Locations_ InSitu::locations() const {
  typedef oops::SampledLocations<ObsTraits> SampledLocations_;
  return Locations_(SampledLocations_(obsdb_.locations()));
}

// -----------------------------------------------------------------------------

void InSitu::print(std::ostream & os) const {
  if (obsdb_.comm().rank() == 0) {
    os << "AQ in situ observation operator";
  }
}

// -----------------------------------------------------------------------------

}  // namespace aq
