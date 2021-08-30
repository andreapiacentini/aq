/*
 * (C) Copyright 2019-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef AQ_MODEL_GETVALUESAQ_H_
#define AQ_MODEL_GETVALUESAQ_H_

#include <memory>
#include <ostream>
#include <string>

#include "eckit/config/LocalConfiguration.h"

#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "oops/aq/AqFortran.h"
#include "oops/aq/LocationsAQ.h"

namespace eckit {
  class Configuration;
}

namespace aq {
  class GomAQ;
  class GeometryAQ;
  class StateAQ;

/// \brief used for getting state values at observation locations
// -----------------------------------------------------------------------------
class GetValuesAQ : public util::Printable,
                    private util::ObjectCounter<GetValuesAQ> {
 public:
  static const std::string classname() {return "aq::GetValuesAQ";}

/// \brief saves all locations \p locs to use during filling GeoVaLs
  GetValuesAQ(const GeometryAQ &, const LocationsAQ & locs, const eckit::Configuration &);
  ~GetValuesAQ() {}

/// \brief fills in \p geovals for all observations in the timeframe (\p t1, \p t2],
/// \p geovals are interpolated trilinearly from \p state at the nearest gridpoints
  void fillGeoVaLs(const StateAQ &, const util::DateTime & t1,
                   const util::DateTime & t2, GomAQ &) const;

 private:
  void print(std::ostream &) const;
  LocationsAQ locs_;
  eckit::LocalConfiguration conf_;
};
// -----------------------------------------------------------------------------

}  // namespace aq

#endif  // AQ_MODEL_GETVALUESAQ_H_
