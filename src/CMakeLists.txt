#==================
# mimir-api library
#==================

set(MIMIR_SRC
  mimir/MimirPriv.cpp
  mimir/StateMachine.cpp
  mimir/algorithm/TestAlgorithm.cpp
  mimir/algorithm/KinematicVessel.cpp
  mimir/algorithm/FishSchool.cpp
  mimir/algorithm/Leadline.cpp
  mimir/algorithm/PursePlannerFormulation.cpp
  mimir/algorithm/PursePlanner.cpp
  mimir/algorithm/AlgorithmFactory.cpp

  )

set(MIMIR_PUBLIC
  ${PROJECT_SOURCE_DIR}/include/mimir/IAlgorithm.hpp)

add_library(mimir-iface INTERFACE)
target_include_directories(mimir-iface INTERFACE
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)
set_target_properties(mimir-iface
  PROPERTIES
  EXPORT_NAME "${PROJECT_NAME}"
  #PUBLIC_HEADER "${MIMIR_PUBLIC}" # Not supported until 3.19, need install files too..
  )

add_library(${PROJECT_NAME}::mimir-iface ALIAS mimir-iface)

# Because PUBLIC_HEADER property on header only target is not supported yet
install(FILES
  ${MIMIR_PUBLIC}
  COMPONENT dev
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/mimir)

install(TARGETS mimir-iface
  EXPORT "${PROJECT_NAME}"
  LIBRARY
  COMPONENT dev
  DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE
  COMPONENT dev
  DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME
  COMPONENT runtime
  DESTINATION ${CMAKE_INSTALL_BINDIR}
  PUBLIC_HEADER
  COMPONENT dev
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/mimir)

add_library(mimir-api OBJECT ${MIMIR_SRC})

target_link_libraries(mimir-api PUBLIC
  ${MIMIR_DDS_TARGETS}
  Boost::boost
  Boost::log
  casadi
  gnuplot-iostream
  mimir-iface
  yaml-cpp
  yavl-cpp)

target_include_directories(mimir-api
PRIVATE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

target_compile_definitions(mimir-api PUBLIC BOOST_LOG_DYN_LINK)

#=======================================
# Executable object
#=======================================

set(PROG_EXEC_FILES
  mimir/program/Config.cpp
  mimir/program/Convenience.cpp
  mimir/program/Options.cpp)

add_library(mimir-exec OBJECT
  ${PROG_EXEC_FILES})

target_include_directories(mimir-exec
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:mimir>)

target_link_libraries(mimir-exec
  PUBLIC
  Boost::boost
  Boost::log
  Boost::program_options
  casadi
  yaml-cpp
  yavl-cpp
  $<$<AND:$<CXX_COMPILER_ID:GNU>,$<VERSION_LESS:$<CXX_COMPILER_VERSION>,9.0>>:stdc++fs> # std::filesystem linking..
  $<$<BOOL:${UNIX}>:pthread> # mingw_stdthreads?
  )
target_compile_definitions(mimir-exec PUBLIC BOOST_LOG_DYN_LINK)


#=======================================
# Executable: mimir

set(HUMAN_APP_NAME mimir)
set(MIMIR_APP_NAME ${HUMAN_APP_NAME} PARENT_SCOPE)
set(APP_NAME mimir-app)

set(PROG_FILES
  programs/mimir.cpp)

add_executable(${APP_NAME} ${PROG_FILES})
target_link_libraries(${APP_NAME}
  PRIVATE
  ${MIMIR_DDS_TARGETS}
  mimir-api
  mimir-exec
  yavl-cpp
  $<$<AND:$<CXX_COMPILER_ID:GNU>,$<VERSION_LESS:$<CXX_COMPILER_VERSION>,9.0>>:stdc++fs> # std::filesystem linking..
  )
target_compile_definitions(${APP_NAME} PUBLIC BOOST_LOG_DYN_LINK)

target_include_directories(${APP_NAME}
  PRIVATE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

set_target_properties(${APP_NAME}
  PROPERTIES
  DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX}
  OUTPUT_NAME mimir
  INSTALL_RPATH "/opt/opensplice-rts/lib;/opt/opensplice-hde/lib;/opt/ospl/lib;/opt/fkin/lib"
  COMPILE_DEFINITIONS MIMIR_INSTALL_PREFIX="${MIMIR_INSTALL_PREFIX}"
  )

install(TARGETS ${APP_NAME}
  RUNTIME
    COMPONENT runtime
    DESTINATION ${CMAKE_INSTALL_BINDIR})

if(MIMIR_CONAN)
  set_target_properties(${APP_NAME}
  PROPERTIES
  INSTALL_RPATH "\$ORIGIN/../lib")
endif()


#=======================================
# Executable: testing gnuplot-iostream

# need boost: iostreams, system, filesystem
# and apt-get install gnuplot gnuplot-iostream-dev

if(MIMIR_WITH_GNUPLOT)
  set(APP_NAME mimir-testing)

  set(PROG_FILES
    programs/testing.cpp)

  add_executable(${APP_NAME} ${PROG_FILES})
  target_link_libraries(${APP_NAME}
    PRIVATE
    ${MIMIR_DDS_TARGETS}
    casadi
    gnuplot-iostream
    mimir-api
    mimir-exec
    yaml-cpp
    yavl-cpp
    $<$<AND:$<CXX_COMPILER_ID:GNU>,$<VERSION_LESS:$<CXX_COMPILER_VERSION>,9.0>>:stdc++fs>
    )
  target_compile_definitions(${APP_NAME} PUBLIC BOOST_LOG_DYN_LINK)

  target_include_directories(${APP_NAME}
    PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

  set_target_properties(${APP_NAME}
    PROPERTIES
    DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX}
    OUTPUT_NAME testing
    INSTALL_RPATH "/opt/opensplice-rts/lib;/opt/opensplice-hde/lib;/opt/ospl/lib;/opt/fkin/lib"
    )

endif()


#========================================
# Executable: checking yavl-cpp

set(APP_NAME yaml-testing)

set(PROG_FILES
  programs/checker.cpp)

add_executable(${APP_NAME} ${PROG_FILES})
target_link_libraries(${APP_NAME}
  PRIVATE
  casadi
  mimir-exec
  yaml-cpp
  yavl-cpp
  )
#target_compile_definitions(${APP_NAME} PUBLIC BOOST_LOG_DYN_LINK)

target_include_directories(${APP_NAME}
  PRIVATE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

set_target_properties(${APP_NAME}
  PROPERTIES
  DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX}
  OUTPUT_NAME checker
  )
