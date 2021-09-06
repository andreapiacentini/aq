/*
 * (C) Copyright 2019-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef AQ_GETVALUES_H_
#define AQ_GETVALUES_H_

#include <memory>
#include <ostream>
#include <string>

#include "eckit/config/LocalConfiguration.h"

#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "aq/Locations.h"

namespace eckit {
  class Configuration;
}

namespace aq {
  class GeoVals;
  class Geometry;
  class State;

/// \brief used for getting state values at observation locations
// -----------------------------------------------------------------------------
class GetValue : public util::Printable,
                    private util::ObjectCounter<GetValue> {
 public:
  static const std::string classname() {return "aq::GetValue";}

/// \brief saves all locations \p locs to use during filling GeoVaLs
  GetValue(const Geometry &, const Locations & locs, const eckit::Configuration &);
  ~GetValue() {}

/// \brief fills in \p geovals for all observations in the timeframe (\p t1, \p t2],
/// \p geovals are interpolated trilinearly from \p state at the nearest gridpoints
  void fillGeoVaLs(const State &, const util::DateTime & t1,
                   const util::DateTime & t2, GeoVals &) const;

 private:
  void print(std::ostream &) const;
  Locations locs_;
  eckit::LocalConfiguration conf_;
};
// -----------------------------------------------------------------------------

}  // namespace aq

#endif  // AQ_GETVALUES_H_
