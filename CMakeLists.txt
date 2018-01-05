cmake_minimum_required(VERSION 2.8)

get_filename_component(FIPS_ROOT_DIR "../fips" ABSOLUTE)
include("${FIPS_ROOT_DIR}/cmake/fips.cmake")


fips_setup(PROJECT oryol-delaunay)
fips_add_subdirectory(src)
fips_finish()
