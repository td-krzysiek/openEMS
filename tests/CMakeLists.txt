# C++ unit tests, run with `ctest` in the build directory.
#
# Each test is a self-contained executable linking the openEMS library that was
# just built. Add a new one to TESTS below, no further wiring needed.

set(TESTS
  test_arraylib
  test_engine_plan
  test_cell_model
  test_cartesian_subgrid
  test_cartesian_subgrid_mesh
  test_cartesian_subgrid_material
  test_cartesian_subgrid_pec
  test_cartesian_subgrid_ports
  test_cartesian_subgrid_boundaries
  test_cartesian_subgrid_mur
  test_cartesian_subgrid_pml
  test_cartesian_subgrid_factor
  test_cartesian_subgrid_bound
  test_cartesian_subgrid_curl
  test_cartesian_subgrid_device
  test_cartesian_subgrid_device_ports
  test_cartesian_subgrid_device_pml
  test_cartesian_subgrid_device_precision
  test_cartesian_subgrid_device_patterns
)

foreach(test ${TESTS})
  add_executable(${test} ${test}.cpp)
  # the tests use the headers of the source tree, not the installed ones
  target_include_directories(${test} PRIVATE ${openEMS_SOURCE_DIR})
  target_link_libraries(${test} openEMS)

  if(test MATCHES "^test_cartesian_subgrid_device" AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|IntelLLVM")
    target_compile_options(${test} PRIVATE -ffp-contract=off)
  endif()
  if(ENABLE_HIP AND NOT OPENEMS_GPU_PLUGINS AND (OPENEMS_HIP_AS_CUDA OR CMAKE_HIP_PLATFORM STREQUAL "nvidia"))
    target_compile_definitions(${test} PRIVATE OPENEMS_SUBGRID_TEST_CUDA=1)
  endif()
  add_test(NAME ${test} COMMAND ${test})
endforeach()

# Previously failing interface DC stress; keep the original 2e-4 field gate.
add_test(NAME test_cartesian_subgrid_device_dc COMMAND test_cartesian_subgrid_device_ports --dc-stress)
