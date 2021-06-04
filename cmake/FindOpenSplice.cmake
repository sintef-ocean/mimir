# Find module for Vortex OpenSplice DDS (OSS edition)
#
# If successful, the script defines the following entities.
#
# IMPORTED targets:
#
#     OpenSplice::isocpp
#     OpenSplice::isocpp2
#
# Variables:
#
#     OpenSplice_FOUND
#     OpenSplice_IDLPP_BINARY
#     OpenSplice_INCLUDE_DIRS
#     OpenSplice_LIBRARIES     (= "OpenSplice::isocpp2")
#     OpenSplice_VERSION
#
# Functions:
#
#     OpenSplice_generate_isocpp(<IDL file> <output directory> <variable for source file names>)
#     OpenSplice_generate_isocpp2(<IDL file> <output directory> <variable for source file names>)
#
# ==============================================================================

# Checks that `osplHomeDir` is an OpenSplice installation directory for
# build configuration `config`, and if so, returns the OpenSplice version
# in a variable whose name is given by `versionVar`.
function(OpenSplice_check_dir osplHomeDir config versionVar)
    if(WIN32)
        set(releaseScript "release.bat")
    else()
        set(releaseScript "release.com")
    endif()
    file(READ "${osplHomeDir}/${releaseScript}" releaseScriptContents)
    string(REGEX MATCH
        "OpenSplice HDE Release ([0-9]+\\.[0-9]+\\.[0-9]+)OSS For ([a-z0-9_]+\\.[a-z0-9_]+(-[a-z0-9_]+)?),"
        match
        "${releaseScriptContents}")
    if("${CMAKE_MATCH_2}" STREQUAL "${config}")
        set(${versionVar} "${CMAKE_MATCH_1}" PARENT_SCOPE)
    else()
        set(${versionVar} "NOTFOUND" PARENT_SCOPE)
    endif()
endfunction()

# Configures all targets for build configuration `config`, assuming that
# `homeDir` actually contains the installed files for this configuration.
function(OpenSplice_configure_targets config homeDir csEnabled)
    string(TOUPPER "${config}" configUpper)
    set(SUFFIX "_${configUpper}")

    if(csEnabled AND EXISTS "${homeDir}/bin/dcpssacsAssembly.dll")
        set(addCSharp TRUE)
    else()
        set(addCSharp FALSE)
    endif()
    if(NOT TARGET "OpenSplice::isocpp2")
        add_library("OpenSplice::ddskernel" SHARED IMPORTED)
        add_library("OpenSplice::ddsi2" SHARED IMPORTED)
        add_library("OpenSplice::spliced" SHARED IMPORTED)
        add_library("OpenSplice::durability" SHARED IMPORTED)

        set_property(TARGET "OpenSplice::ddskernel" PROPERTY INTERFACE_INCLUDE_DIRECTORIES
            "${homeDir}/include"
            "${homeDir}/include/sys")
        add_library("OpenSplice::isocpp" SHARED IMPORTED)
        set_property(TARGET "OpenSplice::isocpp" PROPERTY INTERFACE_INCLUDE_DIRECTORIES
            "${homeDir}/include/dcps/C++/isocpp"
            "${homeDir}/include/dcps/C++/SACPP")
        set_property(TARGET "OpenSplice::isocpp" PROPERTY INTERFACE_LINK_LIBRARIES
            "OpenSplice::ddskernel")
        add_library("OpenSplice::isocpp2" SHARED IMPORTED)
        set_property(TARGET "OpenSplice::isocpp2" PROPERTY INTERFACE_INCLUDE_DIRECTORIES
            "${homeDir}/include/dcps/C++/isocpp2"
            "${homeDir}/include/dcps/C++/SACPP")
        set_property(TARGET "OpenSplice::isocpp2" PROPERTY INTERFACE_LINK_LIBRARIES
            "OpenSplice::ddskernel")
    endif()
    if(addCSharp AND NOT TARGET "OpenSplice::sacs")
        add_library("OpenSplice::sacs" SHARED IMPORTED)
    endif()
    set_property(TARGET "OpenSplice::ddskernel" APPEND PROPERTY IMPORTED_CONFIGURATIONS "${configUpper}")
    set_property(TARGET "OpenSplice::isocpp" APPEND PROPERTY IMPORTED_CONFIGURATIONS "${configUpper}")
    set_property(TARGET "OpenSplice::isocpp2" APPEND PROPERTY IMPORTED_CONFIGURATIONS "${configUpper}")
    set_property(TARGET "OpenSplice::ddsi2" APPEND PROPERTY IMPORTED_CONFIGURATIONS "${configUpper}")
    set_property(TARGET "OpenSplice::spliced" APPEND PROPERTY IMPORTED_CONFIGURATIONS "${configUpper}")
    set_property(TARGET "OpenSplice::durability" APPEND PROPERTY IMPORTED_CONFIGURATIONS "${configUpper}")
    if(WIN32)
        set_target_properties("OpenSplice::ddskernel" PROPERTIES
            IMPORTED_IMPLIB${SUFFIX} "${homeDir}/lib/ddskernel.lib"
            IMPORTED_LOCATION${SUFFIX} "${homeDir}/bin/ddskernel.dll")
          set_target_properties("OpenSplice::ddsi2" PROPERTIES
            IMPORTED_IMPLIB${SUFFIX} "${homeDir}/lib/ddsi2.lib"
            IMPORTED_LOCATION${SUFFIX} "${homeDir}/bin/ddsi2.dll")
          set_target_properties("OpenSplice::spliced" PROPERTIES
            IMPORTED_IMPLIB${SUFFIX} "${homeDir}/lib/spliced.lib"
            IMPORTED_LOCATION${SUFFIX} "${homeDir}/bin/spliced.dll")
          set_target_properties("OpenSplice::durability" PROPERTIES
            IMPORTED_IMPLIB${SUFFIX} "${homeDir}/lib/durability.lib"
            IMPORTED_LOCATION${SUFFIX} "${homeDir}/bin/durability.dll")
        set_target_properties("OpenSplice::isocpp" PROPERTIES
            IMPORTED_IMPLIB${SUFFIX} "${homeDir}/lib/dcpsisocpp.lib"
            IMPORTED_LOCATION${SUFFIX} "${homeDir}/bin/dcpsisocpp.dll")
        set_target_properties("OpenSplice::isocpp2" PROPERTIES
            IMPORTED_IMPLIB${SUFFIX} "${homeDir}/lib/dcpsisocpp2.lib"
            IMPORTED_LOCATION${SUFFIX} "${homeDir}/bin/dcpsisocpp2.dll")
        if(addCSharp)
            set_property(TARGET "OpenSplice::sacs" APPEND PROPERTY IMPORTED_CONFIGURATIONS "${configUpper}")
            set_target_properties("OpenSplice::sacs" PROPERTIES
                IMPORTED_COMMON_LANGUAGE_RUNTIME${SUFFIX} "CSharp"
                IMPORTED_LOCATION${SUFFIX} "${homeDir}/bin/dcpssacsAssembly.dll")
        endif()
    else()
        set_target_properties("OpenSplice::ddskernel" PROPERTIES
            IMPORTED_LOCATION${SUFFIX} "${homeDir}/lib/libddskernel.so"
            IMPORTED_NO_SONAME${SUFFIX} TRUE)
          set_target_properties("OpenSplice::ddsi2" PROPERTIES
            IMPORTED_LOCATION${SUFFIX} "${homeDir}/lib/libddsi2.so"
            IMPORTED_NO_SONAME${SUFFIX} TRUE)
          set_target_properties("OpenSplice::spliced" PROPERTIES
            IMPORTED_LOCATION${SUFFIX} "${homeDir}/lib/libspliced.so"
            IMPORTED_NO_SONAME${SUFFIX} TRUE)
          set_target_properties("OpenSplice::durability" PROPERTIES
            IMPORTED_LOCATION${SUFFIX} "${homeDir}/lib/libdurability.so"
            IMPORTED_NO_SONAME${SUFFIX} TRUE)
        set_target_properties("OpenSplice::isocpp" PROPERTIES
            IMPORTED_LOCATION${SUFFIX} "${homeDir}/lib/libdcpsisocpp.so"
            IMPORTED_NO_SONAME${SUFFIX} TRUE)
        set_target_properties("OpenSplice::isocpp2" PROPERTIES
            IMPORTED_LOCATION${SUFFIX} "${homeDir}/lib/libdcpsisocpp2.so"
            IMPORTED_NO_SONAME${SUFFIX} TRUE)
    endif()
endfunction()


if(NOT DEFINED OSPL_HOME)
    if(DEFINED ENV{OSPL_HOME})
        file(TO_CMAKE_PATH $ENV{OSPL_HOME} OSPL_HOME)
    endif(DEFINED ENV{OSPL_HOME})
endif(NOT DEFINED OSPL_HOME)

if(WIN32)
    get_property(csharp_enabled GLOBAL PROPERTY ENABLED_LANGUAGES)
    if(csharp_enabled MATCHES "CSharp")
        set(csharp_enabled TRUE)
    else()
      set(csharp_enabled FALSE)
    endif()
    if(CMAKE_SIZEOF_VOID_P EQUAL 8 OR (csharp_enabled AND CMAKE_CSharp_COMPILER_ARCHITECTURE_ID STREQUAL "x64"))
        set(OpenSplice_config "x86_64.win64")
    else()
        set(OpenSplice_config "x86.win32")
    endif()
else()
    set(csharp_enabled FALSE)

    if("${CMAKE_CXX_COMPILER}" MATCHES ".*clang.*")
        set(suffx "_clang")
    endif()

    if(DEFINED ANDROID)
      string(PREPEND suffx "_android")
    endif()

    if(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
      set(OpenSplice_config "x86_64.linux")
    elseif(CMAKE_SYSTEM_PROCESSOR STREQUAL "armv7l")
      set(OpenSplice_config "armv7l.linux")
    elseif(CMAKE_SYSTEM_PROCESSOR STREQUAL "armv-v7a") # Android
      set(OpenSplice_config "armv7l.linux")
    elseif(CMAKE_SYSTEM_PROCESSOR STREQUAL "armv7-a") # Android
      set(OpenSplice_config "armv7l.linux")
    elseif(CMAKE_SYSTEM_PROCESSOR STREQUAL "aarch64")
      set(OpenSplice_config "armv8.linux")
    elseif(CMAKE_SYSTEM_PROCESSOR STREQUAL "arm64-v8a") # Android
      set(OpenSplice_config "armv8.linux")
    else()
      set(OpenSplice_config "x86.linux")
    endif()

    string(APPEND OpenSplice_config "${suffx}")

endif(WIN32)

find_path(OpenSplice_HOME_RELEASE
    NAMES
        "include/dcps/C++/isocpp2/dds/dds.hpp"
    HINTS
        "${OSPL_HOME}"
        /opt/opensplice
    PATH_SUFFIXES
        "${OpenSplice_config}"
        "OpenSplice/${OpenSplice_config}"
    NO_CMAKE_FIND_ROOT_PATH
)
mark_as_advanced(OpenSplice_HOME_RELEASE)
find_path(OpenSplice_HOME_DEBUG
    NAMES
        "include/dcps/C++/isocpp2/dds/dds.hpp"
    HINTS
        "${OSPL_HOME}"
        /opt/opensplice
        /opt/opensplice-hde-dbg
    PATH_SUFFIXES
        "${OpenSplice_config}-dev"
        "OpenSplice/${OpenSplice_config}-dev"
    NO_CMAKE_FIND_ROOT_PATH
)
mark_as_advanced(OpenSplice_HOME_DEBUG)
OpenSplice_check_dir("${OpenSplice_HOME_RELEASE}" "${OpenSplice_config}" OpenSplice_versionRelease)
OpenSplice_check_dir("${OpenSplice_HOME_DEBUG}" "${OpenSplice_config}-dev" OpenSplice_versionDebug)

if(OpenSplice_versionRelease OR OpenSplice_versionDebug)
    if(OpenSplice_versionRelease)
        OpenSplice_configure_targets("release" "${OpenSplice_HOME_RELEASE}" ${csharp_enabled})
        if(OpenSplice_versionDebug STREQUAL OpenSplice_versionRelease)
            OpenSplice_configure_targets("debug" "${OpenSplice_HOME_DEBUG}" ${csharp_enabled})
        endif()
        set(OpenSplice_VERSION "${OpenSplice_versionRelease}")
        set(OpenSplice_home "${OpenSplice_HOME_RELEASE}")
    else() # OpenSplice_versionDebug is set
        OpenSplice_configure_targets("debug" "${OpenSplice_HOME_DEBUG}" ${csharp_enabled})
        set(OpenSplice_VERSION "${OpenSplice_versionDebug}")
        set(OpenSplice_home "${OpenSplice_HOME_DEBUG}")
    endif()
    get_property(
        OpenSplice_INCLUDE_DIRS
        TARGET "OpenSplice::isocpp2"
        PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
    set(OpenSplice_LIBRARIES "OpenSplice::isocpp2")
    set(OpenSplice_IDLPP_BINARY "${OpenSplice_home}/bin/idlpp")

    if(WIN32)
        set(OpenSplice_idlpp_wrapper "${CMAKE_BINARY_DIR}/FindOpenSplice_idlpp_wrapper.bat")
        file(WRITE "${OpenSplice_idlpp_wrapper}"
            "${OpenSplice_home}/release.bat && \"${OpenSplice_IDLPP_BINARY}\" %*\n")
    else()
        file(WRITE "${CMAKE_BINARY_DIR}/FindOpenSplice_tmp/FindOpenSplice_idlpp_wrapper"
            "#!/bin/bash\nsource \"${OpenSplice_home}/release.com\" && \"${OpenSplice_IDLPP_BINARY}\" \"$@\"\n")
        file(COPY "${CMAKE_BINARY_DIR}/FindOpenSplice_tmp/FindOpenSplice_idlpp_wrapper"
            DESTINATION "${CMAKE_BINARY_DIR}"
            FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ)
        file(REMOVE_RECURSE "${CMAKE_BINARY_DIR}/FindOpenSplice_tmp")
        set(OpenSplice_idlpp_wrapper "${CMAKE_BINARY_DIR}/FindOpenSplice_idlpp_wrapper")
    endif()

    function(OpenSplice_generate_isocpp idlFile outputDir sourceFilesVar)
        get_filename_component(baseName "${idlFile}" NAME_WE)
        set(sourceFiles
            "${outputDir}/${baseName}.cpp"
            "${outputDir}/${baseName}.h"
            "${outputDir}/${baseName}Dcps.cpp"
            "${outputDir}/${baseName}Dcps.h"
            "${outputDir}/${baseName}Dcps_impl.cpp"
            "${outputDir}/${baseName}Dcps_impl.h"
            "${outputDir}/${baseName}SplDcps.cpp"
            "${outputDir}/${baseName}SplDcps.h"
            "${outputDir}/${baseName}_DCPS.hpp")
        add_custom_command(
            OUTPUT ${sourceFiles}
            COMMAND "${OpenSplice_idlpp_wrapper}" "-S" "-l" "isocpp" "-d" "${outputDir}" "${idlFile}"
            MAIN_DEPENDENCY "${idlFile}"
            VERBATIM)
        set(${sourceFilesVar} ${sourceFiles} PARENT_SCOPE)
    endfunction()

    function(OpenSplice_generate_isocpp2 idlFile outputDir sourceFilesVar)
        get_filename_component(baseName "${idlFile}" NAME_WE)
        set(sourceFiles
            "${outputDir}/${baseName}.cpp"
            "${outputDir}/${baseName}.h"
            "${outputDir}/${baseName}SplDcps.cpp"
            "${outputDir}/${baseName}SplDcps.h"
            "${outputDir}/${baseName}_DCPS.hpp")
        add_custom_command(
            OUTPUT ${sourceFiles}
            COMMAND "${OpenSplice_idlpp_wrapper}" "-S" "-l" "isocpp2" "-d" "${outputDir}" "${idlFile}"
            MAIN_DEPENDENCY "${idlFile}"
            VERBATIM)
        set(${sourceFilesVar} ${sourceFiles} PARENT_SCOPE)
    endfunction()

    if(TARGET "OpenSplice::sacs")
        function(OpenSplice_generate_sacs idlFile outputDir sourceFilesVar)
            get_filename_component(baseName "${idlFile}" NAME_WE)
            set(sourceFiles
                "${outputDir}/${baseName}.cs"
                "${outputDir}/I${baseName}Dcps.cs"
                "${outputDir}/${baseName}Dcps.cs"
                "${outputDir}/${baseName}SplDcps.cs")
            add_custom_command(
                OUTPUT ${sourceFiles}
                COMMAND "${OpenSplice_idlpp_wrapper}" "-S" "-l" "cs" "-d" "${outputDir}" "${idlFile}"
                MAIN_DEPENDENCY "${idlFile}"
                VERBATIM)
            set(${sourceFilesVar} ${sourceFiles} PARENT_SCOPE)
        endfunction()
    endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OpenSplice DEFAULT_MSG
    OpenSplice_IDLPP_BINARY
    OpenSplice_INCLUDE_DIRS
    OpenSplice_LIBRARIES
    OpenSplice_VERSION)
