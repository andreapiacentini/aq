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
