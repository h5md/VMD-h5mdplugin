set(CMAKE_BUILD_TYPE_INIT "Release")

if(DEFINED CMAKE_C_COMPILER_ID)
  if(CMAKE_C_COMPILER_ID MATCHES "GNU|Intel|Clang")
    set(CMAKE_C_FLAGS_INIT "-Wall -std=c99 -pedantic -fPIC")
  endif()
endif()

if(CMAKE_C_PLATFORM_ID STREQUAL "Linux")
  set(CMAKE_SHARED_LINKER_FLAGS_INIT "-Wl,--no-undefined -Wl,--as-needed")
  set(CMAKE_SHARED_LINKER_FLAGS_RELEASE_INIT "-Wl,-s")

  # set appropriate RPATH on installed binaries as well as in build tree,
  # see http://www.vtk.org/Wiki/CMake_RPATH_handling
  #
  # use, i.e. don't skip, the full RPATH for the build tree
  set(CMAKE_SKIP_BUILD_RPATH  FALSE)
  # when building, don't use the install RPATH already
  # (but later on when installing)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
  # add the automatically determined parts of the RPATH
  # which point to directories outside the build tree to the install RPATH
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
endif()
