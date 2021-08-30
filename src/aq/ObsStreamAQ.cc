/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "model/ObsStreamAQ.h"

#include <vector>

#include "eckit/config/Configuration.h"
#include "model/GomAQ.h"
#include "model/ObsBias.h"
#include "model/ObsSpaceAQ.h"
#include "model/ObsVecAQ.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

// -----------------------------------------------------------------------------
namespace aq {
// -----------------------------------------------------------------------------
static ObsOpMaker<ObsStreamAQ> makerStream_("Stream");
// -----------------------------------------------------------------------------

ObsStreamAQ::ObsStreamAQ(const ObsSpaceAQ & odb, const eckit::Configuration & config)
  : obsdb_(odb), varin_(std::vector<std::string>{config.getString("mod var")})
{
  oops::Log::trace() << "ObsStreamAQ created." << std::endl;
}

// -----------------------------------------------------------------------------

void ObsStreamAQ::simulateObs(const GomAQ & gom, ObsVecAQ & ovec,
                              const ObsBias & bias) const {
  aq_stream_equiv_f90(gom.toFortran(), ovec.toFortran(), bias.stream());
}

// -----------------------------------------------------------------------------

std::unique_ptr<LocationsAQ> ObsStreamAQ::locations() const {
  return obsdb_.locations();
}

// -----------------------------------------------------------------------------

void ObsStreamAQ::print(std::ostream & os) const {
  os << "AQ Stream observation operator TL/AD";
}

// -----------------------------------------------------------------------------

}  // namespace aq
