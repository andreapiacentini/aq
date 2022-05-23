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

#include "aq/ObsInSitu.h"

#include <vector>

#include "aq/aq_insitu_interface.h"
#include "aq/GeoVals.h"
#include "aq/ObsAuxControl.h"
#include "aq/ObsSpace.h"
#include "aq/ObsVec.h"

#include "eckit/config/Configuration.h"

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

std::unique_ptr<Locations> InSitu::locations() const {
  return obsdb_.locations();
}

// -----------------------------------------------------------------------------

void InSitu::print(std::ostream & os) const {
  if (obsdb_.comm().rank() == 0) {
    os << "AQ in situ observation operator";
  }
}

// -----------------------------------------------------------------------------

}  // namespace aq
