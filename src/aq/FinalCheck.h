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

#ifndef AQ_FINALCHECK_H_
#define AQ_FINALCHECK_H_

#include <memory>
#include <ostream>

#include "eckit/config/LocalConfiguration.h"

#include "aq/Traits.h"

#include "oops/base/Variables.h"
#include "oops/interface/ObsFilterBase.h"
#include "oops/util/Printable.h"

namespace aq {
  class GeoVals;
  template <typename DATATYPE> class ObsData;
  class ObsDiagnostics;
  class ObsSpace;
  class ObsVec;

class FinalCheck : public oops::interface::ObsFilterBase<ObsTraits> {
 public:
  FinalCheck(const ObsSpace &, const eckit::Configuration &,
            std::shared_ptr<ObsData<int> >, std::shared_ptr<ObsData<float> >): novars_() {}

  void preProcess() override {}
  void priorFilter(const GeoVals &) override {}
  void postFilter(const GeoVals &,
                  const ObsVec &,
                  const ObsVec &,
                  const ObsDiagnostics &) override {}
  void checkFilterData(const oops::FilterStage filterStage) override {}

  oops::Variables requiredVars() const override {return novars_;}
  oops::Variables requiredHdiagnostics() const override {return novars_;}

 private:
  void print(std::ostream &) const override {}
  const oops::Variables novars_;
};

}  // namespace aq

#endif  // AQ_FINALCHECK_H_
