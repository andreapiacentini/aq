/*
 * (C) Copyright 2020-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef AQ_ANALYTICINIT_H_
#define AQ_ANALYTICINIT_H_

#include "eckit/config/LocalConfiguration.h"

namespace aq {
  class Locations;
  class GeoVals;

/// AnalyticInit class fills GeoVaLs with analytic formulae
/// Options: baroclinic instability and large vortices
class AnalyticInit {
 public:
  explicit AnalyticInit(const eckit::Configuration &);
  void fillGeoVaLs(const Locations &, GeoVals &) const;

 private:
  const eckit::LocalConfiguration config_;
};

}  // namespace aq

#endif  // AQ_ANALYTICINIT_H_
