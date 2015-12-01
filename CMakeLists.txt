cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
set(CMAKE_USER_MAKE_RULES_OVERRIDE "cmake/platform.cmake")

project(VMD-h5mdplugin C)

find_package(HDF5 REQUIRED COMPONENTS "C" "HL")
# find_package(VMD REQUIRED)
# FIXME e.g., https://github.com/gromacs/gromacs/blob/master/cmake/FindVMD.cmake
set(GMX_VMD_PLUGIN_PATH "$ENV{VMDDIR}/plugins/include" CACHE STRING "Path to VMD plugin headers")
message(STATUS "VMD plugin path: ${GMX_VMD_PLUGIN_PATH}")

include_directories("${HDF5_INCLUDE_DIR}")
include_directories("${GMX_VMD_PLUGIN_PATH}")

#
# add H5MD library
#
add_library(h5md # SHARED
    libh5md.c
)

#
# add H5MD plugin
#
add_library(h5mdplugin SHARED
    h5mdplugin.c
)
target_link_libraries(h5md
    ${HDF5_LIBRARIES}
    # math library
    m
)

# unset "lib" prefix
set_target_properties(h5mdplugin PROPERTIES PREFIX "")

target_link_libraries(h5mdplugin
    h5md
    ${HDF5_LIBRARIES}
)

install(TARGETS h5mdplugin
  DESTINATION "$ENV{VMDDIR}/plugins/LINUXAMD64/molfile"
)

#
# add test
#
enable_testing()

add_executable(h5mdtest
    h5mdtest.c
    h5mdplugin.c
)
target_link_libraries(h5mdtest
    h5md
    ${HDF5_LIBRARIES}
    # math library
    m
)

add_test(h5md/particles/full_vmd_structure
    h5mdtest h5md "${CMAKE_SOURCE_DIR}/samples/full_vmd_structure.h5"
)
add_test(h5md/particles/half_vmd_structure
    h5mdtest h5md "${CMAKE_SOURCE_DIR}/samples/half_vmd_structure.h5"
)
add_test(h5md/particles/no_vmd_structure
    h5mdtest h5md "${CMAKE_SOURCE_DIR}/samples/no_vmd_structure.h5"
)
