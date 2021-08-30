/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "aq/ObsOpBaseAQ.h"

#include "eckit/config/Configuration.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

namespace aq {

// -----------------------------------------------------------------------------

ObsOpFactory::ObsOpFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    oops::Log::error() << name << " already registered in aq::ObsOpFactory." << std::endl;
    ABORT("Element already registered in ObsOpFactory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

ObsOpBaseAQ * ObsOpFactory::create(const ObsSpaceAQ & odb,
                                   const eckit::Configuration & conf) {
  oops::Log::trace() << "ObsOpBaseAQ::create starting" << std::endl;
  oops::Log::debug() << "ObsOpBaseAQ::create conf" << conf << std::endl;
  const std::string id = conf.getString("obs type");
  typename std::map<std::string, ObsOpFactory*>::iterator
    jloc = getMakers().find(id);
  if (jloc == getMakers().end()) {
    oops::Log::error() << id << " does not exist in observation operator factory." << std::endl;
    ABORT("Element does not exist in ObsOpFactory.");
  }
  ObsOpBaseAQ * ptr = jloc->second->make(odb, conf);
  oops::Log::trace() << "ObsOpBaseAQ::create done" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------

}  // namespace aq
