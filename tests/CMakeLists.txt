# C++ unit tests, run with `ctest` in the build directory.
#
# Each test is a self-contained executable linking the openEMS library that was
# just built. Add a new one to TESTS below, no further wiring needed.

set(TESTS
  test_arraylib
  test_engine_plan
  test_cell_model
)

foreach(test ${TESTS})
  add_executable(${test} ${test}.cpp)
  # the tests use the headers of the source tree, not the installed ones
  target_include_directories(${test} PRIVATE ${openEMS_SOURCE_DIR})
  target_link_libraries(${test} openEMS)
  add_test(NAME ${test} COMMAND ${test})
endforeach()

add_executable(test_observation_profile test_observation_profile.cpp)
target_include_directories(test_observation_profile PRIVATE ${openEMS_SOURCE_DIR})
target_link_libraries(test_observation_profile openEMS)
add_test(NAME observation_profile_on COMMAND test_observation_profile
  ${CMAKE_CURRENT_BINARY_DIR}/observation-profile-test.jsonl on)
add_test(NAME observation_profile_off COMMAND test_observation_profile
  ${CMAKE_CURRENT_BINARY_DIR}/observation-profile-disabled.jsonl off)
