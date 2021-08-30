/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "model/AqTraits.h"
#include "oops/runs/Run.h"
#include "test/interface/GeoVaLs.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  test::GeoVaLs<aq::AqObsTraits> tests;
  return run.execute(tests);
}

