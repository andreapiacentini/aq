/*
 * (C) Copyright 2020-2020 UCAR
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
    // Qui devo avere una cosa del genere
    // obsdb_.getdb("ObsValue", yobs);
    oops::Log::info() << filterConfig << std::endl;
    aq_obsvec_threshold_check_f90(hofx.toFortran(), yobs.toFortran(),
      (*qcflags_).toFortran(), filterConfig);
    size_t inflate = 0;
    size_t ireject = 0;
    // for (size_t jj = 0; jj < yobs.size(); ++jj) {
    //   if (std::abs(yobs[jj] - hofx[jj]) > options_.threshold) {
    //     // inflate obs error variance
    //     if (options_.inflation.value() != boost::none) {
    //       (*obserr_)[jj] *= *options_.inflation.value();
    //       ++inflate;
    //     // or reject observation
    //     } else {
    //       (*qcflags_)[jj] = 1;
    //       ++ireject;
    //     }
    //   }
    // }
    oops::Log::info() << "BackgroundCheck::postFilter rejected = " << ireject
                      << ", inflated = " << inflate << std::endl;
  }
}

// -----------------------------------------------------------------------------
void BackgroundCheck::print(std::ostream & os) const {
  os << "AQ Background check with absolute threshold " << options_.threshold;
}

}  // namespace aq

