# sfh_check_doc_tools()
#
# Internal helper macro to checks if doxygen, python and sphinx are installed
#
macro(sfh_check_doc_tools)

   find_package(Doxygen REQUIRED)
   find_program(SPHINX_EXECUTABLE NAMES sphinx-build DOC "Path to sphinx-build executable")

   set(_doc_requirements
      "Documentation requires doxygen, python3, pip and several other modules.\n\
 To disable documentation: specify '-D ADD_DOC_TARGET=OFF' to cmake")

   if(NOT DOXYGEN_FOUND)
      message( FATAL_ERROR "Doxygen not found. ${_doc_requirements}")
   endif()

   find_program(Python_EXECUTABLE NAMES python python3 DOC "Path to python executable")
   if(Python_EXECUTABLE)
       execute_process(
         COMMAND ${Python_EXECUTABLE} --version
         ERROR_VARIABLE _PYTHON_VERSION
         OUTPUT_VARIABLE _PYTHON_VERSION)
       string(REGEX MATCH "Python 2\..*" _PYTHON_REG "${_PYTHON_VERSION}")
       if(_PYTHON_REG)
         find_program(Python3_EXECUTABLE NAMES python3 DOC "Path to python executable")
         execute_process(
           COMMAND ${Python3_EXECUTABLE} --version
           ERROR_VARIABLE _PYTHON_VERSION
           OUTPUT_VARIABLE _PYTHON_VERSION)
       else()
         set(Python3_EXECUTABLE ${Python_EXECUTABLE})
       endif()
       string(REGEX MATCH "Python 3\..*" _PYTHON_REG "${_PYTHON_VERSION}")
       if(NOT _PYTHON_REG)
         message(FATAL_ERROR "Wrong python version: ${_PYTHON_VERSION}")
       else()
         set(Python3_FOUND 1)# PARENT_SCOPE)
       endif()
   endif()

   if(Python3_FOUND)
      message(STATUS "Found python: ${Python3_EXECUTABLE}")
   else()
      message(FATAL_ERROR "Python3 interpreter not found. ${_doc_requirements}")
   endif()

   if(SPHINX_EXECUTABLE)
      message(STATUS "Found sphinx: ${SPHINX_EXECUTABLE}")
   else()
      message(FATAL_ERROR "Sphinx executable not found. ${_doc_requirements}")
   endif()

endmacro()


# sfh_check_pip()
#
# Checks whether pip is installed
#
function(sfh_check_pip)

   sfh_check_doc_tools()

   # Check if pip is installed
   execute_process(
      COMMAND ${Python3_EXECUTABLE} -c "import pip"
      RESULT_VARIABLE PIP_EXIT_CODE
      OUTPUT_QUIET
      ERROR_QUIET)

   if(NOT PIP_EXIT_CODE EQUAL 0)
      message(FATAL_ERROR "Python3 with pip not found.")
   endif()
endfunction()


# sfh_doc_python_requirements(_requirements)
#
# Checks if contents pip ${_requirements} file modules are installed
#
function(sfh_doc_python_requirements _requirements)

   sfh_check_pip()

   # Check python modules
   execute_process(
      COMMAND ${Python3_EXECUTABLE} -m pip freeze -r ${_requirements}
      RESULT_VARIABLE _req_status
      ERROR_VARIABLE _req_freeze_warning
      OUTPUT_QUIET
      )

   # pip freeze complained, checking required modules one by one
   if(_req_freeze_warning)
      message(WARNING "pip freeze warning: ${_req_freeze_warning}")
      file(STRINGS requirements.txt _python_requirements)
      foreach(_py_req IN LISTS _python_requirements)
         string(REGEX MATCH "^[A-Za-z0-9]+[A-Za-z0-9_-]*[A-za-z0-9]*" _py_mod ${_py_req})
         string(LENGTH "${_py_mod}" _py_mod_length)
         if(_py_mod_length EQUAL 0)
            continue()
         endif()
         message(STATUS "Checking python module: ${_py_req}")
         execute_process(
            COMMAND ${Python3_EXECUTABLE} -m pip show ${_py_mod}
            RESULT_VARIABLE _req_status
            OUTPUT_QUIET
            ERROR_QUIET
            )
         if(NOT _req_status EQUAL 0)
            list(APPEND _py_req_missing ${_py_req})
         endif()
      endforeach()
      list(LENGTH _py_req_missing _py_modules_missing)
      if(_py_modules_missing GREATER 0)
         string(REPLACE ";" " " _py_req_missing "${_py_req_missing}")
         set(_py_install_req "\nInstall modules with 'python3 -m pip install -r ${CMAKE_CURRENT_SOURCE_DIR}/requirements.txt'")
         message(FATAL_ERROR "Missing python3 modules: ${_py_req_missing}. ${_py_install_req}")
      endif()
   else()
      message(STATUS "Found required python modules: OK")
   endif()

endfunction()
