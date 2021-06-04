# The documentation is built with several targets. It makes use of doxygen to generate
# xml, which is used by sphinx, breathe, exhale to generate html documentation.
# Exhale is used to autogenerate api doc
# sphinx doc source directory in build_dir/sphinx/source to contain conf.py and .rst/.md
# The documentation creation also uses emacs to convert from org-mode to .rst

include(SfhDocHelpers.cmake)

sfh_check_doc_tools()
sfh_doc_python_requirements("${CMAKE_SOURCE_DIR}/docs/requirements.txt")
find_program(DOC2DASH_EXECUTABLE doc2dash DOC "Path to doc2dash executable")
find_program(EMACS_EXECUTABLE emacs "Path to emacs executable")

# --------------------------------------------------------------------
# Doxygen xml generation :: common settings

set(DOXYGEN_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doxygen)

# Global generation settings
set(DOXYGEN_GENERATE_HTML NO)
set(DOXYGEN_HTML_OUTPUT html)
set(DOXYGEN_USE_MATHJAX YES) # To disable latex usage in html.
set(DOXYGEN_GENERATE_XML YES)
set(DOXYGEN_GENERATE_LATEX NO)
if(WIN32 AND NOT UNIX)
  set(DOXYGEN_PLANTUML_JAR_PATH C:/ProgramData/chocolatey/lib/plantuml/tools/plantuml.jar)
else()
  set(DOXYGEN_PLANTUML_JAR_PATH /usr/share/plantuml/plantuml.jar) # optional
endif()
message(STATUS "Plantuml: ${DOXYGEN_PLANTUML_JAR_PATH}")

# Suggested settings by Exhale docs
set(DOXYGEN_EXTRACT_ALL YES)
set(DOXYGEN_FULL_PATH_NAMES YES) # Are stripped by exhale
set(DOXYGEN_ENABLE_PREPROCESSING YES)
set(DOXYGEN_MACRO_EXPANSION YES)
set(DOXYGEN_EXPAND_ONLY_PREDEF NO)
set(DOXYGEN_SKIP_FUNCTION_MACROS NO)
set(DOXYGEN_XML_PROGRAMLISTING YES)

# Enable restructuredText in doxygen
list(APPEND DOXYGEN_ALIASES "rst=\\verbatim embed:rst:")
list(APPEND DOXYGEN_ALIASES "endrst=\\endverbatim")
list(APPEND DOXYGEN_ALIASES "rststar=\\verbatim embed:rst:leading-asterisk")
list(APPEND DOXYGEN_ALIASES "endrststar=\\endverbatim")
list(APPEND DOXYGEN_ALIASES "rstinline=\\verbatim embed:rst:inline ")
list(APPEND DOXYGEN_ALIASES "endrstinline=\\endverbatim")

# Tweak parsing to avoid some issues with Sphinx generation
list(APPEND DOXYGEN_PREDEFINED DOXYGEN_SHOULD_SKIP_THIS)
list(APPEND DOXYGEN_PREDEFINED "sequence=sequence_") # IDL: avoids Sphinx compilation error

# More doxygen graphs does not affect XML
# Do not extract all, like cpp files and more
# Tweak what you want included in doxygen html
# Note: It seem to only affect non-XML export
set(DOXYGEN_HAVE_DOT YES)
set(DOXYGEN_UML_LOOK YES)
set(DOXYGEN_CLASS_DIAGRAMS YES)
set(DOXYGEN_CLASS_GRAPH YES)
set(DOXYGEN_COLLABORATION_GRAPH YES)
set(DOXYGEN_HIDE_UNDOC_MEMBERS NO)
set(DOXYGEN_HIDE_UNDOC_CLASSES NO)
set(DOXYGEN_HIDE_FRIEND_COMPOUNDS NO)
set(DOXYGEN_HIDE_IN_BODY_DOCS YES)


# Clean doxygen folder
#
add_custom_target(doxy-clean-xml
  COMMAND ${CMAKE_COMMAND} -E remove_directory ${DOXYGEN_OUTPUT_DIRECTORY}
  COMMENT "DELETE generated doxygen"
  )

# --------------------------------------------------------------------
# Doxygen generation for public api
# --------------------------------------------------------------------

list(APPEND DOXYGEN_EXCLUDE_PATTERNS "*.cpp") # Do not parse cpp files, avoids some issues

set(DOXYGEN_XML_OUTPUT xml)
set(DOXYGEN_HTML_OUTPUT html)
set(DOXYGEN_STRIP_FROM_PATH ${CMAKE_SOURCE_DIR}/include)
set(MY_SOURCES "../include")

doxygen_add_docs(doxygen "${MY_SOURCES}"
  COMMENT "Generating docs with Doxygen")
# CMake>=3.16: USE_STAMP_FILE to avoid building all the time?

add_dependencies(doxygen doxy-clean-xml)

# --------------------------------------------------------------------
# Doxygen generation for src, including idl qml
# --------------------------------------------------------------------

set(DOXYGEN_XML_OUTPUT xml_priv)
set(DOXYGEN_HTML_OUTPUT html_priv)
set(DOXYGEN_STRIP_FROM_PATH ${CMAKE_SOURCE_DIR})

#list(APPEND DOXYGEN_FILE_PATTERNS "*.idl") # Comment to disable idl structs in docs
list(APPEND DOXYGEN_FILE_PATTERNS "*.hpp") # Comment to not parse hpp files from src

set(MY_SOURCES "../src/mimir")
list(APPEND MY_SOURCES "../src")

doxygen_add_docs(doxygen-priv "${MY_SOURCES}"
  COMMENT "Generating docs with Doxygen")

add_dependencies(doxygen-priv doxygen)

# --------------------------------------------------------------------
# Sphinx html generation
# --------------------------------------------------------------------

set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/sphinx)
set(SPHINX_SOURCE ${SPHINX_BUILD}/source)
set(SPHINX_INDEX_FILE ${SPHINX_BUILD}/html/index.html)
set(SPHINX_LATEX_FILE ${SPHINX_BUILD}/latex/${PROJECT_NAME}.tex)
set(SPHINX_DOCSET_FILE ${SPHINX_BUILD}/${PROJECT_NAME}.docset/Contents/Info.plist)
set(SPHINX_MAN_FILE ${SPHINX_BUILD}/man/testproject.1)

string(TIMESTAMP THE_YEAR '%Y')

if (WITH_API_DOC)
  set(MIMIR_WITH_API_DOC "'t") # or 'nil if to disable, latter is discouraged
else()
  set(MIMIR_WITH_API_DOC "'nil")
  if (NOT WIN32)
    set(MIMIR_LATEX ON)
  else()
    message(WARNING "Latex documentation output is not yet supported for windows")
  endif()
endif()


configure_file(conf.py.in ${SPHINX_SOURCE}/conf.py @ONLY)

add_custom_target(doc-org DEPENDS ${SPHINX_SOURCE}/index.rst)
add_dependencies(doc-org doxygen-priv) # doxygen doxygen-priv)

file(GLOB_RECURSE THE_ORG_FILES
  LIST_DIRECTORIES false
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  CONFIGURE_DEPENDS
  org/*.org)

# Convert org files to restructuredText using emacs
add_custom_command(
  OUTPUT ${SPHINX_SOURCE}/index.rst always.outdated
  COMMAND ${CMAKE_COMMAND} -E make_directory ${SPHINX_SOURCE}/rst
  COMMAND ${EMACS_EXECUTABLE} -Q --batch
  -L ${CMAKE_CURRENT_SOURCE_DIR}
  --eval "(setq out-dir \"${SPHINX_SOURCE}\")"
  --eval "(setq with-api-doc ${MIMIR_WITH_API_DOC})"
  -l emacs-load.el
  -f jemacs-export-org-doc
  index.org
  COMMAND ${EMACS_EXECUTABLE} -Q --batch
  -L ${CMAKE_CURRENT_SOURCE_DIR}
  --eval "(setq out-dir \"${SPHINX_SOURCE}/rst\")"
  -l emacs-load.el
  -f jemacs-export-org-doc
  ${THE_ORG_FILES}
  COMMAND ${CMAKE_COMMAND} -E copy_directory static ${SPHINX_SOURCE}/static
  MAIN_DEPENDENCY index.org
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  COMMENT "Exporting org files to restructuredText"
  VERBATIM
  )


# Prepare sphinx root source directory
# Clean api build every time
# Copy static resources
# Copy markdown files
add_custom_target(doc-sources
  COMMAND ${CMAKE_COMMAND} -E remove_directory ${SPHINX_SOURCE}/api
  COMMAND ${CMAKE_COMMAND} -E remove_directory ${SPHINX_SOURCE}/api_private
  COMMAND ${CMAKE_COMMAND} -E remove_directory ${SPHINX_SOURCE}/org
  COMMAND ${CMAKE_COMMAND} -E remove_directory ${SPHINX_BUILD}/html
  COMMAND ${CMAKE_COMMAND} -E copy_directory static ${SPHINX_SOURCE}/static
  COMMAND ${CMAKE_COMMAND} -E copy_directory _templates ${SPHINX_SOURCE}/_templates
  COMMAND ${CMAKE_COMMAND} -E copy_directory _ext ${SPHINX_SOURCE}/_ext
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/data/images/FHF-black-txt.png ${SPHINX_SOURCE}/static
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/data/icons/logo.png ${SPHINX_SOURCE}/static
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/data/icons/icon.ico ${SPHINX_SOURCE}/static
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  COMMENT "DELETE generated api docs and COPY doc sources to sphinx root dir"
  )

add_custom_target(doc-html DEPENDS ${SPHINX_INDEX_FILE})
add_dependencies(doc-html doc-org doc-sources)# doxygen doxygen-priv)

# Generate the sphinx html
add_custom_command(
  OUTPUT ${SPHINX_INDEX_FILE} always.outdated
  COMMAND ${SPHINX_EXECUTABLE} -b html
  -Dbreathe_projects.${PROJECT_NAME}=${DOXYGEN_OUTPUT_DIRECTORY}/xml
  -Dbreathe_projects.${PROJECT_NAME}-private=${DOXYGEN_OUTPUT_DIRECTORY}/xml_priv
  ${SPHINX_SOURCE} ${SPHINX_BUILD}/html
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  MAIN_DEPENDENCY ${SPHINX_SOURCE}/conf.py
  COMMENT "Generate docs with Sphinx"
  VERBATIM
  )

install(DIRECTORY ${SPHINX_BUILD}/html
  COMPONENT doc
  DESTINATION ${CMAKE_INSTALL_DOCDIR}
  PATTERN ".doctrees*" EXCLUDE
  PATTERN "_sources*" EXCLUDE
  PATTERN ".buildinfo" EXCLUDE
  )

# Add custom command here to enable build of pdf
# only output latex if api doc is disabled
if(MIMIR_LATEX)

  find_program(LATEXMK_EXECUTABLE latexmk "Path to latexmk executable")
  if(NOT LATEXMK_EXECUTABLE)
    message(WARNING "Latexmk not was found")
    set(LATEX_CMD echo)
  else()
    set(LATEX_CMD ${LATEXMK_EXECUTABLE})
  endif()

  set(MIMIR_LATEX_TARGET doc-latex)
  add_custom_target(doc-latex DEPENDS ${SPHINX_LATEX_FILE})
  add_dependencies(doc-latex doc-org doc-sources)

  set(MIMIR_LATEX_PATCH ${SPHINX_BUILD}/latex/patch.sh)
  set(MIMIR_LATEX_BUILD ${SPHINX_BUILD}/latex/build.sh)
  # write file
  file(WRITE ${MIMIR_LATEX_PATCH} "#!/bin/bash\n")
  file(APPEND ${MIMIR_LATEX_PATCH} "pushd ${SPHINX_BUILD}/latex\n")
  file(APPEND ${MIMIR_LATEX_PATCH} "sed -i 's/\\\\begin{split}\\\\begin{equation}/%/g' ${PROJECT_NAME}.tex\n")
  file(APPEND ${MIMIR_LATEX_PATCH} "sed -i 's/\\\\end{equation}\\\\end{split}/%/g' ${PROJECT_NAME}.tex\n")
  file(APPEND ${MIMIR_LATEX_PATCH} "sed -i 's|equation:rst/\\(.*\\):prob|rst/\\1:prob|g' ${PROJECT_NAME}.tex\n")
  file(APPEND ${MIMIR_LATEX_PATCH} "sed -i 's|equation:rst/\\(.*\\):crit|rst/\\1:crit|g' ${PROJECT_NAME}.tex\n")
  file(APPEND ${MIMIR_LATEX_PATCH} "sed -i 's|equation:rst/\\(.*\\):ass|rst/\\1:ass|g' ${PROJECT_NAME}.tex\n")
  file(APPEND ${MIMIR_LATEX_PATCH} "sed -i 's|equation:rst/\\(.*\\):def|rst/\\1:def|g' ${PROJECT_NAME}.tex\n")
  file(APPEND ${MIMIR_LATEX_PATCH} "popd\n")

  file(WRITE ${MIMIR_LATEX_BUILD} "#!/bin/bash\n")
  file(APPEND ${MIMIR_LATEX_BUILD} "pushd ${SPHINX_BUILD}/latex\n")
  file(APPEND ${MIMIR_LATEX_BUILD} "${LATEX_CMD} -pdf ${PROJECT_NAME}.tex\n")
  file(APPEND ${MIMIR_LATEX_BUILD} "popd\n")

  add_custom_command(
    OUTPUT ${SPHINX_LATEX_FILE} always.outdated
    COMMAND ${SPHINX_EXECUTABLE} -b latex ${SPHINX_SOURCE} ${SPHINX_BUILD}/latex
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/docs/tex/defs.tex ${SPHINX_BUILD}/latex
    COMMAND chmod +x ${MIMIR_LATEX_PATCH}
    COMMAND ${MIMIR_LATEX_PATCH}
    COMMAND chmod +x ${MIMIR_LATEX_BUILD}
    COMMAND ${MIMIR_LATEX_BUILD}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    MAIN_DEPENDENCY ${SPHINX_SOURCE}/conf.py
    COMMENT "Generate latex docs with Sphinx"
    VERBATIM
    )
endif()


# --------------------------------------------------------------------
# Docset generation based on html doc
# --------------------------------------------------------------------

add_custom_target(docset DEPENDS ${SPHINX_DOCSET_FILE})
add_dependencies(docset doc-html)

# Create docset of html
add_custom_command(
  OUTPUT ${SPHINX_DOCSET_FILE} always.outdated
  COMMAND ${DOC2DASH_EXECUTABLE} -j -f -n ${PROJECT_NAME}
  -d ${SPHINX_BUILD} ${SPHINX_BUILD}/html
  COMMAND
    ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/data/icons/icon@1x.png
    ${SPHINX_BUILD}/${PROJECT_NAME}.docset/icon.png
  COMMAND
    ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/data/icons/icon@2x.png
    ${SPHINX_BUILD}/${PROJECT_NAME}.docset/icon@2x.png
  COMMAND
    ${CMAKE_COMMAND} -E tar czf
    ${PROJECT_NAME}_${PROJECT_VERSION}_docset.tar.gz
    ${PROJECT_NAME}.docset
  WORKING_DIRECTORY ${SPHINX_BUILD}
  COMMENT "Generate docset with doc2dash"
  VERBATIM
  )

set(FEED_HOST http://localhost:20000)
configure_file(docset_feed.xml.in ${SPHINX_BUILD}/${PROJECT_NAME}.xml)
# If you host the folder with the docset.tar.gz and xml it will be available
# Quick and dirty: python -m http.server 20000

# Docset target is not to be installed

# --------------------------------------------------------------------
# Meta target for all docs
# --------------------------------------------------------------------

# If you add ALL to custom target, doc will always build when you call cmake --build .
add_custom_target(doc)
add_dependencies(doc doc-html ${MIMIR_LATEX_TARGET} docset)

add_custom_target(docs)
add_dependencies(docs doc)
