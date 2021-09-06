/*
 * (C) Copyright 2019-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <memory>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "aq/GetValues.h"

#include "oops/util/Logger.h"

#include "aq/Geometry.h"
#include "aq/GeoVals.h"
#include "aq/Locations.h"
#include "aq/State.h"

namespace aq {

// -----------------------------------------------------------------------------
/// Constructor, destructor
// -----------------------------------------------------------------------------
GetValues::GetValues(const Geometry &, const Locations & locs,
                         const eckit::Configuration & conf)
    : locs_(locs), conf_(conf)
{
  oops::Log::trace() << "GetValues constructor with config " << conf_ << std::endl;
}

// -----------------------------------------------------------------------------
/// Get state values at observation locations
// -----------------------------------------------------------------------------
void GetValues::fillGeoVaLs(const State & state, const util::DateTime & t1,
                              const util::DateTime & t2, GeoVals & geovals) const
{
  oops::Log::trace() << "GetValues::fillGeoVaLs start" << std::endl;
  // the below call is an example if one wanted a different interpolation type
  const std::string interpType = conf_.getString("interpolation type", "default");

  if (interpType == "default" || (interpType.compare(0, 8, "default_") == 0)) {
    oops::Log::trace() << state.fields() << std::endl;
    aq_getvalues_interp_f90(locs_, state.fields().toFortran(), t1, t2, geovals.toFortran());
  } else {
    std::string err_message("interpolation type option " + interpType + " not supported");
    throw eckit::BadValue(err_message, Here());
  }
  oops::Log::trace() << "GetValues::fillGeoVaLs done" << std::endl;
}

// -----------------------------------------------------------------------------
void GetValues::print(std::ostream & os) const {
  os << "AQ GetValues";
}
// -----------------------------------------------------------------------------

}  // namespace aq
