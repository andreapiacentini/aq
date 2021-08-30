/*
 * (C) Copyright 2018-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "aq/AqTraits.h"
// #include "oops/aq/instantiateAqChangeVarFactory.h"
#include "oops/runs/HofX3D.h"
#include "oops/runs/Run.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  // aq::instantiateAqChangeVarFactory();
  oops::HofX3D<aq::AqTraits, aq::AqObsTraits> hofx;
  return run.execute(hofx);
}
