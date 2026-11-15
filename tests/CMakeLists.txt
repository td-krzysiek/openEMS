# C++ unit tests, run with `ctest` in the build directory.
#
# Each test is a self-contained executable linking the openEMS library that was
# just built. Add a new one to TESTS below, no further wiring needed.

set(TESTS
  test_arraylib
  test_engine_plan
  test_cell_model
  test_setup_profile
  test_lossymetal_setup_profile
  test_lossymetal_extent_cache
)

foreach(test ${TESTS})
  add_executable(${test} ${test}.cpp)
  # the tests use the headers of the source tree, not the installed ones
  target_include_directories(${test} PRIVATE ${openEMS_SOURCE_DIR})
  target_link_libraries(${test} openEMS)
  add_test(NAME ${test} COMMAND ${test})
endforeach()

# Manual integration diagnostic used by check_lossymetal_setup.py.
add_executable(inspect_lossymetal_setup inspect_lossymetal_setup.cpp)
target_include_directories(inspect_lossymetal_setup PRIVATE ${openEMS_SOURCE_DIR})
target_link_libraries(inspect_lossymetal_setup openEMS)
