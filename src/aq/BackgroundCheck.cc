/*
 * (C) Copyright 2020-2020 UCAR
 * (C) Copyright 2021-2022 CERFACS.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "aq/BackgroundCheck.h"
#include "aq/aq_obsvec_interface.h"
#include "aq/ObsVec.h"

#include "aq/Traits.h"

// -----------------------------------------------------------------------------
namespace aq {
// -----------------------------------------------------------------------------
static oops::interface::FilterMaker<ObsTraits, BackgroundCheck> makerBackgroundCheck_(
    "Background Check");

// -----------------------------------------------------------------------------
BackgroundCheck::BackgroundCheck(const ObsSpace & obsdb, const Parameters_ & parameters,
           std::shared_ptr<ObsData<int> > qcflags, std::shared_ptr<ObsData<float> > obserr)
  : obsdb_(obsdb), options_(parameters), qcflags_(qcflags), obserr_(obserr), novars_()
{
}

// -----------------------------------------------------------------------------
void BackgroundCheck::postFilter(const GeoVals &,
                                 const ObsVec & hofx,
                                 const ObsVec &,
                                 const ObsDiagnostics &) {
  if (obsdb_.comm().rank() == 0) {
    ObsVec yobs(obsdb_, "ObsValue");
    oops::Log::info() << "BackgroundCheck::postFilter operations " << std::endl;
    eckit::LocalConfiguration filterConfig(options_.toConfiguration());
    oops::Log::info() << filterConfig << std::endl;
    aq_obsvec_threshold_check_f90(hofx.toFortran(), yobs.toFortran(),
      (*qcflags_).toFortran(), filterConfig);
  }
}

// -----------------------------------------------------------------------------
void BackgroundCheck::print(std::ostream & os) const {
  os << "AQ Background check with absolute threshold " << options_.threshold;
}

}  // namespace aq

