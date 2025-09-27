pybind11_add_module(PyFANS micro.cpp micro.hpp)
target_link_libraries(PyFANS PRIVATE FANS::FANS)

# --- External dependencies ---
target_link_libraries(PyFANS PRIVATE
  Eigen3::Eigen
  nlohmann_json::nlohmann_json
  MPI::MPI_CXX
)
if (TARGET HDF5::HDF5)
  target_link_libraries(PyFANS PRIVATE HDF5::HDF5)
else()
  target_include_directories(PyFANS PRIVATE ${HDF5_INCLUDE_DIRS})
  target_link_libraries(PyFANS PRIVATE ${HDF5_CXX_LIBRARIES})
  target_compile_definitions(PyFANS PRIVATE ${HDF5_DEFINITIONS})
  if (HDF5_IS_PARALLEL)
    target_compile_definitions(PyFANS PRIVATE H5_HAVE_PARALLEL)
  endif()
endif()

if (TARGET FFTW3::FFTW3 AND TARGET FFTW3::MPI)
  target_link_libraries(PyFANS PRIVATE FFTW3::FFTW3 FFTW3::MPI)
else()
  target_include_directories(PyFANS PRIVATE ${FFTW3_INCLUDE_DIRS})
  target_link_libraries(PyFANS PRIVATE ${FFTW3_LIBRARIES})
  target_compile_definitions(PyFANS PRIVATE ${FFTW3_DEFINITIONS})
endif()

add_custom_command(
  TARGET PyFANS
  POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E create_symlink
          $<TARGET_FILE:PyFANS>
          ${CMAKE_CURRENT_SOURCE_DIR}/../test/test_pyfans/$<TARGET_FILE_NAME:PyFANS>
  COMMENT "Create a symlink for FANS python bindings to ${CMAKE_CURRENT_SOURCE_DIR}/../test/"
)
