cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
set(CMAKE_USER_MAKE_RULES_OVERRIDE "cmake/platform.cmake")

project(VMD-h5mdplugin C)

find_package(HDF5 REQUIRED COMPONENTS "C" "HL")

find_program(VMD_EXECUTABLE NAMES vmd PATH_SUFFIXES bin DOC "VMD command")
if(VMD_EXECUTABLE)
    file(STRINGS "${VMD_EXECUTABLE}" VMDDIR REGEX "^[^#].*defaultvmddir=.*$")
    string(REGEX REPLACE "(^.*=\"?|\"$)" "" VMDDIR "${VMDDIR}")
    set(VMD_PLUGIN_INCLUDE_DIR "${VMDDIR}/plugins/include" CACHE STRING "Path to VMD plugin headers")
    set(VMD_PLUGIN_PATH "${VMDDIR}/plugins/LINUXAMD64" CACHE STRING "Path to VMD plugin headers")
    message(STATUS "VMD executable: ${VMD_EXECUTABLE}")
    message(STATUS "VMD plugin path: ${VMD_PLUGIN_PATH}")
    message(STATUS "VMD plugin include dir: ${VMD_PLUGIN_INCLUDE_DIR}")
else()
    message(FATAL_ERROR "VMD not found, please specify with -DVMD_EXECUTABLE")
endif()

include_directories("${HDF5_INCLUDE_DIR}")
include_directories("${VMD_PLUGIN_INCLUDE_DIR}")

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
    ${HDF5_C_LIBRARIES}
    ${HDF5_HL_LIBRARIES}
    # math library
    m
)

SET (CMAKE_SHARED_LINKER_FLAGS
    "${CMAKE_EXE_LINKER_FLAGS} -Wl,-rpath,'$ORIGIN'" )

# unset "lib" prefix
set_target_properties(h5mdplugin PROPERTIES PREFIX "")

target_link_libraries(h5mdplugin
    h5md
    ${HDF5_C_LIBRARIES}
    ${HDF5_HL_LIBRARIES}
)

install(TARGETS h5md h5mdplugin
  DESTINATION "${VMD_PLUGIN_PATH}/molfile"
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
