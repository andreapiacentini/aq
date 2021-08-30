/*
 * (C) Copyright 2019-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef AQ_GETVALUESTLAD_H_
#define AQ_GETVALUESTLAD_H_

#include <memory>
#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"

#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "aq/AqFortran.h"
#include "aq/LocationsAQ.h"

namespace aq {
  class GeometryAQ;
  class GomAQ;
  class IncrementAQ;
  class StateAQ;

/// \brief used for getting state values at observation locations
//  and applying its TL & AD
class GetValuesTLAD : public util::Printable,
                      private util::ObjectCounter<GetValuesTLAD> {
 public:
  static const std::string classname() {return "aq::GetValuesTLAD";}

/// \brief saves all locations \p locs to use during filling GeoVaLs
  GetValuesTLAD(const GeometryAQ &, const LocationsAQ & locs,
                const eckit::Configuration &);
  ~GetValuesTLAD();

  /// \brief fills in \p geovals for all observations in the timeframe (\p t1, \p t2],
  //  \p geovals are interpolated trilinearly from \p state at the nearest gridpoints
  void setTrajectory(const StateAQ & state, const util::DateTime & t1,
                     const util::DateTime & t2, GomAQ & geovals);
  /// \brief fills in \p geovals for all observations in the timeframe (\p t1, \p t2],
  //  as a tangent-linear operator applied to \p inc
  void fillGeoVaLsTL(const IncrementAQ & inc, const util::DateTime & t1,
                     const util::DateTime & t2, GomAQ & geovals) const;
  /// \brief fills in \p inc as adjoint operator applied to \p geovals for all
  //  observations in the timeframe (\p t1, \p t2]
  void fillGeoVaLsAD(IncrementAQ & inc, const util::DateTime & t1,
                     const util::DateTime & t2, const GomAQ & geovals) const;

/// Data
 private:
  void print(std::ostream &) const;
  LocationsAQ locs_;
  F90hmat hmat_;
};
// -----------------------------------------------------------------------------

}  // namespace aq

#endif  // AQ_GETVALUESTLAD_H_
