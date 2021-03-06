/*
 * (C) Copyright 2020-2020 UCAR
 * (C) Copyright 2021-2022 CERFACS.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef AQ_ANALYTICINIT_H_
#define AQ_ANALYTICINIT_H_

#include "oops/interface/AnalyticInitBase.h"

#include "aq/Traits.h"

namespace aq {
  class Locations;
  class GeoVals;

/// Parameters for AQ Analytic init (empty except for analytic init method defined
/// in the base class)
class AnalyticInitParameters : public oops::AnalyticInitParametersBase {
  OOPS_CONCRETE_PARAMETERS(AnalyticInitParameters, AnalyticInitParametersBase)
};

/// AnalyticInit class fills GeoVaLs with analytic formulae
/// Options: baroclinic instability and large vortices
class AnalyticInit : public oops::interface::AnalyticInitBase<ObsTraits> {
 public:
  typedef  AnalyticInitParameters Parameters_;

  explicit AnalyticInit(const Parameters_ &);
  void fillGeoVaLs(const Locations &, GeoVals &) const;

 private:
  Parameters_ options_;
};

}  // namespace aq

#endif  // AQ_ANALYTICINIT_H_
