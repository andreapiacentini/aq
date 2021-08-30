/*
 * (C) Crown Copyright 2021, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef AQ_MODEL_FINALCHECK_H_
#define AQ_MODEL_FINALCHECK_H_

#include <memory>
#include <ostream>

#include "eckit/config/LocalConfiguration.h"

#include "model/AqTraits.h"

#include "oops/base/Variables.h"
#include "oops/interface/ObsFilterBase.h"
#include "oops/util/Printable.h"

namespace aq {
  class GomAQ;
  template <typename DATATYPE> class ObsDataAQ;
  class ObsDiagsAQ;
  class ObsSpaceAQ;
  class ObsVecAQ;

class FinalCheck : public oops::interface::ObsFilterBase<AqObsTraits> {
 public:
  FinalCheck(const ObsSpaceAQ &, const eckit::Configuration &,
             std::shared_ptr<ObsDataAQ<int> >, std::shared_ptr<ObsDataAQ<float> >): novars_() {}

  void preProcess() override {}
  void priorFilter(const GomAQ &) override {}
  void postFilter(const ObsVecAQ &, const ObsDiagsAQ &) override {}

  oops::Variables requiredVars() const override {return novars_;}
  oops::Variables requiredHdiagnostics() const override {return novars_;}

 private:
  void print(std::ostream &) const override {}
  const oops::Variables novars_;
};

}  // namespace aq

#endif  // AQ_MODEL_FINALCHECK_H_
