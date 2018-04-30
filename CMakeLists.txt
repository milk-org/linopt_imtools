project(lib_linopt_imtools_project)

include_directories ("${PROJECT_SOURCE_DIR}/src")
include_directories ("${PROJECT_SOURCE_DIR}/..")
add_library(linopt_imtools SHARED linopt_imtools.c)

install(TARGETS linopt_imtools DESTINATION lib)
install(FILES linopt_imtools.h DESTINATION include)

