/*
 * (C) Copyright 2019-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef AQ_LINEARGETVALUES_H_
#define AQ_LINEARGETVALUES_H_

#include <memory>
#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"

#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "aq/Locations.h"

namespace aq {
  class Geometry;
  class GeoVals;
  class Increment;
  class State;

/// \brief used for getting state values at observation locations
//  and applying its TL & AD
class LinearGetValues : public util::Printable,
                      private util::ObjectCounter<LinearGetValues> {
 public:
  static const std::string classname() {return "aq::LinearGetValues";}

/// \brief saves all locations \p locs to use during filling GeoVaLs
  LinearGetValues(const Geometry &, const Locations & locs,
                const eckit::Configuration &);
  ~LinearGetValues();

  /// \brief fills in \p geovals for all observations in the timeframe (\p t1, \p t2],
  //  \p geovals are interpolated trilinearly from \p state at the nearest gridpoints
  void setTrajectory(const State & state, const util::DateTime & t1,
                     const util::DateTime & t2, GeoVals & geovals);
  /// \brief fills in \p geovals for all observations in the timeframe (\p t1, \p t2],
  //  as a tangent-linear operator applied to \p inc
  void fillGeoVaLsTL(const Increment & inc, const util::DateTime & t1,
                     const util::DateTime & t2, GeoVals & geovals) const;
  /// \brief fills in \p inc as adjoint operator applied to \p geovals for all
  //  observations in the timeframe (\p t1, \p t2]
  void fillGeoVaLsAD(Increment & inc, const util::DateTime & t1,
                     const util::DateTime & t2, const GeoVals & geovals) const;

/// Data
 private:
  void print(std::ostream &) const;
  Locations locs_;
  F90hmat hmat_;
};
// -----------------------------------------------------------------------------

}  // namespace aq

#endif  // AQ_LINEARGETVALUES_H_
