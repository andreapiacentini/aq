/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "model/ObsStreamTLAD.h"

#include <vector>

#include "eckit/config/Configuration.h"
#include "model/AqFortran.h"
#include "model/GomAQ.h"
#include "model/ObsBias.h"
#include "model/ObsBiasIncrement.h"
#include "model/ObsSpaceAQ.h"
#include "model/ObsVecAQ.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

// -----------------------------------------------------------------------------
namespace aq {
// -----------------------------------------------------------------------------
static ObsOpTLADMaker<ObsStreamTLAD> makerStreamTL_("Stream");
// -----------------------------------------------------------------------------

ObsStreamTLAD::ObsStreamTLAD(const ObsSpaceAQ &, const eckit::Configuration & config)
  : varin_(std::vector<std::string>{"x"})
{
  oops::Log::trace() << "ObsStreamTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsStreamTLAD::setTrajectory(const GomAQ &, const ObsBias &) {}

// -----------------------------------------------------------------------------

void ObsStreamTLAD::simulateObsTL(const GomAQ & gom, ObsVecAQ & ovec,
                                  const ObsBiasIncrement & bias) const {
  aq_stream_equiv_tl_f90(gom.toFortran(), ovec.toFortran(), bias.stream());
}

// -----------------------------------------------------------------------------

void ObsStreamTLAD::simulateObsAD(GomAQ & gom, const ObsVecAQ & ovec,
                                  ObsBiasIncrement & bias) const {
  aq_stream_equiv_ad_f90(gom.toFortran(), ovec.toFortran(), bias.stream());
}

// -----------------------------------------------------------------------------

void ObsStreamTLAD::print(std::ostream & os) const {
  os << "AQ Stream observation operator TL/AD";
}

// -----------------------------------------------------------------------------

}  // namespace aq
