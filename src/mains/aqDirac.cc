/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "model/AqTraits.h"
// AQ #include "model/instantiateAqChangeVarFactory.h"
#include "model/instantiateAqLocalizationFactory.h"
#include "oops/runs/Dirac.h"
#include "oops/runs/Run.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  // AQ aq::instantiateAqChangeVarFactory();
  aq::instantiateAqLocalizationFactory();
  oops::Dirac<aq::AqTraits> dir;
  return run.execute(dir);
}
