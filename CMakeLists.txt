cmake_minimum_required(VERSION 3.8)

set(THE_PROJECT_AUTHOR "Joakim Haugen")
set(THE_PROJECT_MAINTAINER "SINTEF Ocean <joakim.haugen@sintef.no>")
set(THE_PROJECT_VENDOR "SINTEF Ocean")

file(STRINGS "${CMAKE_SOURCE_DIR}/version.txt" projectVersion)

find_package(kluster-scripts CONFIG QUIET)

if(kluster-scripts_FOUND)
  # Use git tag and synchronize version.txt
  kluster_get_version_from_git(projectTagVersion ${projectVersion})
  kluster_use_version(${projectTagVersion} "${CMAKE_SOURCE_DIR}/version.txt" projectVersion)
endif()

project(
  mimir
  VERSION ${projectVersion}
  LANGUAGES CXX
  DESCRIPTION "Purse seine deployment decision support")

set(CMAKE_CXX_EXTENSIONS OFF) # Disable compiler flag extension (use -std=c++11 instead of -std=gnu++11)
set(CMAKE_CXX_STANDARD 17) # The default CXX_STANDARD (resolves to a flag -std=c++17)..
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_DEBUG_POSTFIX "_d" CACHE STRING "Generate debug library name with a postfix.")
set(CMAKE_POSITION_INDEPENDENT_CODE ON) # -fPIC by default

option(MIMIR_WITH_GNUPLOT "Enable gnuplot iostream interface" OFF)
option(MIMIR_SKIP_UNOFFICIAL_DEPS "Skip unofficial dependencies in .deb package" OFF)
option(WITH_DOC "Enable building documentation" OFF)
option(WITH_API_DOC "Build API documentation" ON)

# Set some standard variables
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  add_compile_options("-Wall" "-Wextra" "-Wpedantic" "-Werror")
  add_compile_options("-Wno-parentheses")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  add_compile_options("/W4" "/WX")
  add_compile_options("/wd4251")
  add_compile_options("/wd4127")
  add_compile_options("/wd4458")
  add_compile_options("/wd4100")
  add_definitions("-D_SCL_SECURE_NO_WARNINGS" "-D_CRT_SECURE_NO_WARNINGS")
  add_definitions("-D_SILENCE_CXX17_ITERATOR_BASE_CLASS_DEPRECATION_WARNING")
endif()

if(EXISTS ${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
  include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
  conan_basic_setup()
  include(${CMAKE_BINARY_DIR}/conan_paths.cmake)
  set(MIMIR_CONAN ON)
else()
  message(STATUS "Configuring without conan")
  # Use bundled FindOpenSplice
  list(INSERT CMAKE_MODULE_PATH 0 "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
  if(GENERATOR_IS_MULTI_CONFIG)
    message(STATUS "NOTE: Adjusting output directories: Multi-configuration (Release and Debug) builds will replace each other")
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_RELEASE ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_DEBUG ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
  endif()
endif()

if(MIMIR_CONAN)
  set(YAML_FIND MODULE)
  set(MIMIR_INSTALL_PREFIX /opt/sintef)
else()
  set(YAML_FIND CONFIG)
  set(MIMIR_INSTALL_PREFIX /usr/local)
endif()

if(WIN32)
  set(MIMIR_INSTALL_PREFIX "C:/Program Files/${CMAKE_PROJECT_NAME}-${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}")

endif()

find_package(Boost 1.67.0 REQUIRED COMPONENTS log program_options)
find_package(OpenSplice 6.9 REQUIRED)
find_package(casadi 3.5 CONFIG REQUIRED) # Use bundled find package, proper flags and stuff
find_package(yaml-cpp ${YAML_FIND} 0.6 REQUIRED)
find_package(RatatoskIDL CONFIG REQUIRED)


if(MIMIR_WITH_GNUPLOT)
  set(MIMIR_GNUPLOT_BOOST_COMPONENTS iostreams system filesystem)
  find_package(Boost 1.67.0 REQUIRED COMPONENTS ${MIMIR_GNUPLOT_BOOST_COMPONENTS})
  if(NOT MIMIR_CONAN)
    include(CheckIncludeFileCXX)
    CHECK_INCLUDE_FILE_CXX("gnuplot-iostream.h" MIMIR_GNUPLOT_IOSTREAM)
    message(STATUS "GNUPLOT ENABLED, iostream test result: ${MIMIR_GNUPLOT_IOSTREAM}")
  else()
    find_package(gnuplot-iostream MODULE REQUIRED)
    set(MIMIR_GNUPLOT_TARGET gnuplot-iostream::gnuplot-iostream)
    set(MIMIR_GNUPLOT_IOSTREAM ON)
  endif()

  if(MIMIR_GNUPLOT_IOSTREAM)
    message(STATUS "Enabling gnuplot iostream interface")
    add_library(gnuplot-iostream INTERFACE)
    target_link_libraries(gnuplot-iostream
      INTERFACE
      ${MIMIR_GNUPLOT_TARGET}
      Boost::iostreams
      Boost::system
      Boost::filesystem)
    target_compile_definitions(gnuplot-iostream INTERFACE MIMIR_WITH_GNUPLOT=1)
    set(MIMIR_GNUPLOT_DEBIAN_PKGS "gnuplot, libgnuplot-iostream-dev, ")
  else()
    message(FATAL_ERROR "MIMIR_WITH_GNUPLOT requested, but gnuplot-iostream.h not found.")
  endif()
else()
  add_library(gnuplot-iostream INTERFACE)
endif()


# --------------------------------------------------------------------
#  DDS code generation

ratatosk_add_idl_library(mimir-fkin FKIN)
ratatosk_add_idl_library(mimir-ratatosk Ratatosk)

set(MIMIR_DDS_TARGETS "mimir-fkin;mimir-ratatosk")

# --------------------------------------------------------------------
#  Libraries, export targets and install files

# Compile external project yavl-cpp

add_library(yavl-cpp OBJECT
  external/yavl-cpp/src/yavl-cpp/yavl.cpp)

target_compile_features(yavl-cpp INTERFACE cxx_std_17)

target_include_directories(yavl-cpp
  PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/external/yavl-cpp/src/yavl-cpp>
  PUBLIC $<INSTALL_INTERFACE:include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/external/yavl-cpp/include>
)

target_link_libraries(yavl-cpp yaml-cpp)


include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/${CMAKE_PROJECT_NAME})

add_subdirectory(src)

if(WITH_DOC)
  add_subdirectory(docs)
  set(MIMIR_DOC doc)
endif()

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/${PROJECT_NAME}ConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMinorVersion)

configure_file(cmake/projectConfig.cmake.in
  cmake/${PROJECT_NAME}Config.cmake.in @ONLY)

configure_package_config_file(
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/${PROJECT_NAME}Config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/${PROJECT_NAME}Config.cmake
  PATH_VARS INSTALL_CONFIGDIR
  INSTALL_DESTINATION ${INSTALL_CONFIGDIR})

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/${PROJECT_NAME}Config.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/${PROJECT_NAME}ConfigVersion.cmake
  COMPONENT dev
  DESTINATION ${INSTALL_CONFIGDIR})

install(EXPORT ${PROJECT_NAME}
  FILE ${PROJECT_NAME}Targets.cmake
  NAMESPACE "${PROJECT_NAME}::"
  COMPONENT dev
  DESTINATION ${INSTALL_CONFIGDIR})

install(FILES LICENSE
  COMPONENT dev
  DESTINATION ${CMAKE_INSTALL_DATADIR}/${PROJECT_NAME}/dev)

install(FILES LICENSE
  COMPONENT runtime
  DESTINATION ${CMAKE_INSTALL_DATADIR}/${PROJECT_NAME})


# --------------------------------------------------------------------
# Packaging

string(TOLOWER ${CMAKE_CXX_COMPILER_ID} COMPILER)
set(COMPILER_VERSION ${CMAKE_CXX_COMPILER_VERSION})
string(REPLACE "." ";" comp_list "${COMPILER_VERSION}")

list(GET comp_list 0 COMP_MAJOR)
#list(GET comp_list 1 COMP_MINOR) # no need for gnu > 6, but what about other compilers?

set(COMPILER ${COMPILER}-${COMP_MAJOR})

if(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" OR CMAKE_SYSTEM_PROCESSOR STREQUAL "AMD64")
  set(ARCHITECTURE amd64)
elseif(CMAKE_SYSTEM_PROCESSOR STREQUAL "arm64")
  set(ARCHITECTURE arm64)
elseif(CMAKE_SYSTEM_PROCESSOR STREQUAL "armv7l")
  set(ARCHITECTURE armhf)
else()
  set(ARCHITECTURE i386)
endif()


set(THE_PACKAGE_LINUX_INSTALL_PREFIX ${MIMIR_INSTALL_PREFIX})
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "${PROJECT_DESCRIPTION}")
set(CPACK_PACKAGE_VENDOR "${THE_PROJECT_VENDOR}")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README.org")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")

set(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${PROJECT_VERSION_PATCH})
set(CPACK_DEBIAN_PACKAGE_VERSION
  ${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH})

set(CPACK_PACKAGE_FILE_NAME
  ${CMAKE_PROJECT_NAME}-${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}-${ARCHITECTURE})
set(THE_PACKAGE_FILE_NAME ${CPACK_PACKAGE_FILE_NAME}) # needed, since include CPack changes things

set(CPACK_DEB_COMPONENT_INSTALL ON)
set(CPACK_COMPONENTS_ALL runtime dev ${MIMIR_DOC})
set(CPACK_COMPONENT_DEV_DEPENDS runtime ${MIMIR_DOC})

set(CPACK_COMPONENT_DEV_DESCRIPTION "${PROJECT_DESCRIPTION} - development headers")
set(CPACK_COMPONENT_RUNTIME_DESCRIPTION "${PROJECT_DESCRIPTION} - runtime")
if(WITH_DOC)
  set(CPACK_COMPONENT_DOC_DESCRIPTION "${PROJECT_DESCRIPTION} - documentation")
  set(CPACK_DEBIAN_DOC_PACKAGE_DEPENDS "libgomp1")
endif()

if(MIMIR_CONAN)
  list(APPEND CPACK_COMPONENTS_ALL bundle)
  list(APPEND CPACK_COMPONENT_DEV_DEPENDS bundle)
  list(APPEND CPACK_COMPONENT_RUNTIME_DEPENDS bundle)
  set(CPACK_COMPONENT_BUNDLE_DESCRIPTION "${PROJECT_DESCRIPTION} - dependency bundle")
endif()

# -------------------------
# --- WINDOWS packaging ---
if(WIN32 AND NOT UNIX)

  # https://blog.nathanosman.com/2017/11/24/using-windeployqt-with-cpack.html
  set(CMAKE_INSTALL_UCRT_LIBRARIES TRUE)
  include(InstallRequiredSystemLibraries)

  set(CPACK_GENERATOR "NSIS;TGZ")
  set(CPACK_PACKAGE_INSTALL_DIRECTORY ${CMAKE_PROJECT_NAME}-${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR})
  set(CPACK_SOURCE_PACKAGE_FILE_NAME ${CMAKE_PROJECT_NAME}-${VERSION_NUMBER})
  set(CPACK_NSIS_DISPLAY_NAME  "${CPACK_PACKAGE_INSTALL_DIRECTORY}")
  set(CPACK_NSIS_MUI_ICON "${CMAKE_SOURCE_DIR}/data/icons/icon.ico")
  set(CPACK_NSIS_MODIFY_PATH ON)
  set(CPACK_NSIS_CONTACT "${THE_PROJECT_MAINTAINER}")
  set(CPACK_PACKAGE_EXECUTABLES "${MIMIR_APP_NAME}" "${CMAKE_PROJECT_NAME}" )
  set(CPACK_CREATE_DESKTOP_LINKS "${MIMIR_APP_NAME}")
  set(CPACK_NSIS_INSTALLED_ICON_NAME bin\\\\${MIMIR_APP_NAME}.exe)

  set(CPACK_COMPONENT_BUNDLE_HIDDEN ON)
  set(CPACK_COMPONENT_BUNDLE_REQUIRED ON)
  set(CPACK_COMPONENT_RUNTIME_REQUIRED ON)

  if(WITH_DOC)
    set(CPACK_NSIS_MENU_LINKS
      "share/doc/${CMAKE_PROJECT_NAME}/html/index.html" "Mimir Documentation")
  endif()

# ---------------------------
# --- LINUX packaging deb ---
else(WIN32 AND NOT UNIX)

  set(CPACK_GENERATOR "DEB;TGZ")
  set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE ${ARCHITECTURE})
  set(CPACK_DEBIAN_PACKAGE_MAINTAINER "${THE_PROJECT_MAINTAINER}")
  set(CPACK_DEBIAN_PACKAGE_PRIORITY optional)
  set(CPACK_PACKAGING_INSTALL_PREFIX "${THE_PACKAGE_LINUX_INSTALL_PREFIX}")
  set(CPACK_DEBIAN_ENABLE_COMPONENT_DEPENDS ON)

  if(MIMIR_OMP)
    message(STATUS "Debian package: Using OpenMP runtime: ${MIMIR_OMP}")
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    message(STATUS "Debian package: Setting OpenMP runtime to libgomp1, override with MIMIR_OMP")
    set(MIMIR_OMP libgomp1)
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    message(STATUS "Debian package: Setting OpenMP runtime to libomp5, override with MIMIR_OMP")
    set(MIMIR_OMP libomp5)
  endif()

  if(MIMIR_FORTRAN_RT)
    message(STATUS "Debian package: Using Fortran runtime: ${MIMIR_FORTRAN_RT}")
  elseif(COMPILER_CXX_COMPILER_ID STREQUAL "GNU")
    message(STATUS "Debian package: Assuming Fortran runtime: libgfortran, override with MIMIR_FORTRAN_RT")
    if(COMPILER_CXX_COMPILER_VERSION VERSIONLESS 7)
      set(MIMIR_FORTRAN_RT libgfortran3)
    elseif(COMPILER_CXX_COMPILER_VERSION VERSIONEQUAL 7)
      set(MIMIR_FORTRAN_RT libgfortran4)
    else()
      set(MIMIR_FORTRAN_RT libgfortran5)
    endif()
  else()
    message(STATUS "Debian package: Assuming libgfortran5 runtime, override with MIMIR_FORTRAN_RT")
    set(MIMIR_FORTRAN_RT libgfortran5)
  endif()

  if(NOT MIMIR_SKIP_UNOFFICIAL_DEPS)
    set(MIMIR_EXTRA_DEV_DEPS "casadi-dep (= 1.1.0), casadi (>= 3.5), opensplice-hde (>= 6.9.190925)")
    set(MIMIR_EXTRA_RUNTIME_DEPS "casadi-dep (= 1.1.0), casadi (>= 3.5), opensplice-rts (>= 6.9.190925)")
  endif()

  set(CPACK_DEBIAN_DEV_PACKAGE_DEPENDS "\
  ${MIMIR_OMP}, \
  ${MIMIR_FORTRAN_RT}, \
  ${MIMIR_GNUPLOT_DEBIAN_PKGS} \
  libboost-log-dev (>= 1.67.0), \
  libboost-regex-dev (>= 1.67.0), \
  libboost-program-options-dev (>= 1.67.0), \
  libyaml-cpp-dev (>= 0.6), \
  ${MIMIR_EXTRA_DEV_DEPS}")

  if (CMAKE_VERSION VERSION_LESS 3.14)
    set(Boost_VERSION_STRING "${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION}")
  endif()

  set(CPACK_DEBIAN_RUNTIME_PACKAGE_DEPENDS "\
  ${MIMIR_OMP}, \
  ${MIMIR_FORTRAN_RT},
  ${MIMIR_GNUPLOT_DEBIAN_PKGS} \
  libboost-log${Boost_VERSION_STRING}, \
  libboost-program-options${Boost_VERSION_STRING}, \
  libboost-regex${Boost_VERSION_STRING}, \
  libyaml-cpp${yaml-cpp_VERSION_MAJOR}.${yaml-cpp_VERSION_MINOR}, \
  ${MIMIR_EXTRA_RUNTIME_DEPS}")

  # Manpage file generation
  # parse usage.org into mimir.1 using pandoc

  # Man page custom target
  find_program(PANDOC_EXECUTABLE pandoc)
  add_custom_target(manpage DEPENDS ${CMAKE_BINARY_DIR}/mimir.1)

  if(NOT PANDOC_EXECUTABLE)
    set(MIMIR_MAN_PARSE ${CMAKE_COMMAND} -E copy
      ${CMAKE_SOURCE_DIR}/docs/mimir.1 ${CMAKE_BINARY_DIR}/mimir.1)
    set(MIMIR_MAN_PARSE_IN_SOURCE echo "  Warning: manpage only copied: pandoc missing")
    message(STATUS "pandoc not found - mimir.1 manpage may not be in sync with docs/org/usage.org")
  else()
    set(MIMIR_MAN_PARSE
      ${PANDOC_EXECUTABLE} -s ${CMAKE_SOURCE_DIR}/docs/org/usage.org
      -f org -t man -o ${CMAKE_BINARY_DIR}/mimir.1)
    set(MIMIR_MAN_PARSE_IN_SOURCE
      ${PANDOC_EXECUTABLE} -s ${CMAKE_SOURCE_DIR}/docs/org/usage.org
      -f org -t man -o ${CMAKE_SOURCE_DIR}/docs/mimir.1)
  endif()

  add_custom_command(
    DEPENDS ${CMAKE_SOURCE_DIR}/docs/org/usage.org
    OUTPUT ${CMAKE_BINARY_DIR}/mimir.1
    OUTPUT ${CMAKE_BINARY_DIR}/mimir.1.gz
    COMMAND ${MIMIR_MAN_PARSE_IN_SOURCE}
    COMMAND ${MIMIR_MAN_PARSE}
    COMMAND ${CMAKE_COMMAND} -E tar -fcz mimir.1.gz ${CMAKE_SOURCE_DIR}/docs/mimir.1
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Convert usage.org to manpage"
    VERBATIM)


  # $ man /opt/sintef/share/man/mimir
  install(FILES ${CMAKE_BINARY_DIR}/mimir.1.gz
    COMPONENT runtime
    DESTINATION ${CMAKE_INSTALL_MANDIR}/man1)

  # Bundle dependencies directly
  if(MIMIR_CONAN)
    # Dependencies are bundled in the bundle component
    unset(CPACK_DEBIAN_DEV_PACKAGE_DEPENDS)
    unset(CPACK_DEBIAN_RUNTIME_PACKAGE_DEPENDS)
    set(CPACK_DEBIAN_BUNDLE_PACKAGE_DEPENDS "${MIMIR_OMP}, ${MIMIR_FORTRAN_RT}")
  endif()

endif(WIN32 AND NOT UNIX)

# Bundle dependencies in bundle folder
if(MIMIR_CONAN)
  # Install the bundle directory
  install(DIRECTORY
    ${CMAKE_BINARY_DIR}/bundle/
    COMPONENT bundle
    DESTINATION .
    USE_SOURCE_PERMISSIONS)
endif()

include(CPack)

#---------------------------------------------------------------------
# If kluster-scripts found, add upload, deliver, and deploy targets

if(kluster-scripts_FOUND AND UNIX)
  kluster_add_artifactory_target(${MIMIR_DOC})
  kluster_add_aptly_targets(${MIMIR_DOC})
elseif(kluster-scripts_FOUND AND WIN32)
  kluster_add_artifactory_target(${MIMIR_DOC})
else()
  add_custom_target(package_it
    COMMAND cmake --build . --target package --config Release
    COMMENT "Package project ${CMAKE_PROJECT_NAME}"
    VERBATIM)
  if(WITH_DOC)
    add_dependencies(package_it doc)
  endif()
endif()

if(UNIX)
  add_dependencies(package_it manpage)
endif()

if(MIMIR_CONAN AND UNIX)
  find_program(PATCHELF_EXECUTABLE patchelf)

  if(NOT PATCHELF_EXECUTABLE)
    message(WARNING "patchelf tool was not found. Bundled libraries will not have changed runpath")
  else()
    set(MY_ORIGIN "\$ORIGIN")
    file(GLOB MIMIR_BUNDLES
      LIST_DIRECTORIES false
      RELATIVE ${CMAKE_BINARY_DIR}/bundle/lib "${CMAKE_BINARY_DIR}/bundle/lib/*")
    add_custom_target(patch_bundle
      COMMAND ${PATCHELF_EXECUTABLE} --set-rpath ${MY_ORIGIN} ${MIMIR_BUNDLES}
      COMMENT "Patching runpath for bundled libraries"
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bundle/lib
      VERBATIM
      )
    add_dependencies(package_it patch_bundle)
  endif()

endif()
