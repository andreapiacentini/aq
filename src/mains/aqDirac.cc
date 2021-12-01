/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "aq/instantiateAqChangeVarFactory.h"
#include "aq/Traits.h"
#include "oops/runs/Dirac.h"
#include "oops/runs/Run.h"
#include "saber/oops/instantiateCovarFactory.h"
#include "saber/oops/instantiateLocalizationFactory.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  aq::instantiateAqChangeVarFactory();
  saber::instantiateCovarFactory<aq::Traits>();
  saber::instantiateLocalizationFactory<aq::Traits>();
  oops::Dirac<aq::Traits> dir;
  return run.execute(dir);
}
