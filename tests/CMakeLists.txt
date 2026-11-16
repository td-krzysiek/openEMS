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
  test_cartesian_subgrid_device_geometry
  test_cartesian_subgrid_regular_curl
  test_cartesian_subgrid_edge_index
  test_cartesian_subgrid_device
  test_cartesian_subgrid_device_ports
  test_cartesian_subgrid_device_pml
  test_cartesian_subgrid_device_precision
  test_cartesian_subgrid_device_patterns
  test_cartesian_subgrid_device_large
  test_cartesian_subgrid_device_domains
  test_cartesian_subgrid_device_domain_ports
  test_cartesian_subgrid_device_domain_loss
  test_cartesian_subgrid_device_advance
  test_cartesian_subgrid_device_waveforms
  test_cartesian_subgrid_conductivity
  test_cartesian_subgrid_conductivity_guards
  test_cartesian_subgrid_native_conductivity
  test_cartesian_subgrid_conductivity_pml
  test_cartesian_subgrid_surface
  test_cartesian_subgrid_surface_guards
  test_cartesian_subgrid_integral_roundoff
  test_cartesian_subgrid_opaque
  test_cartesian_subgrid_opaque_faces
  test_cartesian_subgrid_native_surface
  test_cartesian_subgrid_shaped_surface
  test_cartesian_subgrid_conductor_union
  test_cartesian_subgrid_native_path_ports
)

foreach(test ${TESTS})
  add_executable(${test} ${test}.cpp)
  # the tests use the headers of the source tree, not the installed ones
  target_include_directories(${test} PRIVATE ${openEMS_SOURCE_DIR})
  target_link_libraries(${test} openEMS)

  if((test MATCHES "^test_cartesian_subgrid_device" OR test STREQUAL "test_cartesian_subgrid_regular_curl") AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|IntelLLVM")
    target_compile_options(${test} PRIVATE -ffp-contract=off)
  endif()
  if(ENABLE_HIP AND NOT OPENEMS_GPU_PLUGINS AND (OPENEMS_HIP_AS_CUDA OR CMAKE_HIP_PLATFORM STREQUAL "nvidia"))
    target_compile_definitions(${test} PRIVATE OPENEMS_SUBGRID_TEST_CUDA=1)
  endif()
  add_test(NAME ${test} COMMAND ${test})
  if(test STREQUAL "test_cartesian_subgrid_native_surface" OR test STREQUAL "test_cartesian_subgrid_conductor_union" OR test STREQUAL "test_cartesian_subgrid_shaped_surface")
    # Volume LossyMetal is optional in the independently versioned CSXCAD API.
    set_tests_properties(${test} PROPERTIES SKIP_RETURN_CODE 77)
  endif()
endforeach()

# Previously failing interface DC stress; keep the original 2e-4 field gate.
add_test(NAME test_cartesian_subgrid_device_dc COMMAND test_cartesian_subgrid_device_ports --dc-stress)
