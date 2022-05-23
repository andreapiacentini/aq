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
