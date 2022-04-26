/*
 * (C) Copyright 2020-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef AQ_BACKGROUNDCHECK_H_
#define AQ_BACKGROUNDCHECK_H_

#include <memory>
#include <ostream>

#include "aq/Traits.h"

#include "oops/base/Variables.h"
#include "oops/generic/ObsFilterParametersBase.h"
#include "oops/interface/ObsFilterBase.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace aq {
  class GeoVals;
  template <typename DATATYPE> class ObsData;
  class ObsSpace;
  class ObsDiagnostics;
  class ObsVec;

/// Parameters for AQ BackgroundCheck
/// background check: all obs for which {|y-H(x)| < threshold} pass QC
class BackgroundCheckParameters : public oops::ObsFilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(BackgroundCheckParameters, ObsFilterParametersBase)

 public:
  /// threshold for background check
  oops::RequiredParameter<double> threshold{"threshold", this};

  /// optional inflation factor: if this parameter is present, obs error stddev
  /// for obs that don't pass the check is multiplied by the specified factor.
  /// Otherwise, obs that don't pass the check are rejected.
  oops::OptionalParameter<double> inflation{"inflate obs error", this};
};

/// Simple background check: all obs for which {|y-H(x)| < threshold} pass QC
class BackgroundCheck : public oops::interface::ObsFilterBase<ObsTraits> {
 public:
  typedef BackgroundCheckParameters Parameters_;

  BackgroundCheck(const ObsSpace &, const Parameters_ &,
                  std::shared_ptr<ObsData<int> >, std::shared_ptr<ObsData<float> >);

  void preProcess() override {}
  void priorFilter(const GeoVals &) override {}
  void postFilter(const GeoVals &,
                  const ObsVec &,
                  const ObsVec &,
                  const ObsDiagnostics &) override;
  void checkFilterData(const oops::FilterStage filterStage) override {}

  oops::Variables requiredVars() const override {return novars_;}
  oops::Variables requiredHdiagnostics() const override {return novars_;}

 private:
  void print(std::ostream & os) const override;

  const ObsSpace & obsdb_;
  Parameters_ options_;
  std::shared_ptr<ObsData<int> > qcflags_;   // QC flags
  std::shared_ptr<ObsData<float> > obserr_;  // obs error stddev
  const oops::Variables novars_;
};

}  // namespace aq

#endif  // AQ_BACKGROUNDCHECK_H_
