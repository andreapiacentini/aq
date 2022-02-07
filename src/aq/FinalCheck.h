/*
 * (C) Crown Copyright 2021, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
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
