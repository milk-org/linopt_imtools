# library name
set(LIBNAME "milklinoptimtools")
set(SRCNAME "linopt_imtools")

message("")
message(" SRCNAME = ${SRCNAME} -> LIBNAME = ${LIBNAME}")

set(SOURCEFILES
	${SRCNAME}.c)

set(INCLUDEFILES
	${SRCNAME}.h)


# DEFAULT SETTINGS
# Do not change unless needed
# =====================================================================

project(lib_${LIBNAME}_project)

include_directories ("${PROJECT_SOURCE_DIR}/src")
include_directories ("${PROJECT_SOURCE_DIR}/..")


# Library can be compiled from multiple source files
# Convention: the main souce file is named <libname>.c
#
add_library(${LIBNAME} SHARED ${SOURCEFILES})
target_include_directories(${LIBNAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)
target_link_libraries(${LIBNAME} PRIVATE CLIcore)

if(USE_CUDA)
find_package(CUDA 9.0 REQUIRED)
message(STATUS "Found CUDA ${CUDA_VERSION_STRING} at ${CUDA_TOOLKIT_ROOT_DIR}")
message("---- CUDA_INCLUDE_DIRS   = ${CUDA_INCLUDE_DIRS}")
message("---- CUDA_LIBRARIES      = ${CUDA_LIBRARIES} ${CUDA_cusolver_LIBRARY}" )
message("---- CUDA_CUBLAS_LIBRARIES    = ${CUDA_CUBLAS_LIBRARIES}" )
target_include_directories(${LIBNAME} PUBLIC ${CUDA_INCLUDE_DIRS})
target_link_libraries(${LIBNAME} PRIVATE ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES} ${CUDA_cusolver_LIBRARY})
target_compile_options(${LIBNAME} PUBLIC -DHAVE_CUDA)
endif(USE_CUDA)



if(USE_MAGMA)
find_package(PkgConfig REQUIRED)
pkg_check_modules(MAGMA REQUIRED magma)
message("---- MAGMA_LIBRARY_DIRS =  ${MAGMA_LIBRARY_DIRS}")
message("---- MAGMA_LIBRARIES    =  ${MAGMA_LIBRARIES}" )
message("---- MAGMA_CFLAGS_OTHER =  ${MAGMA_CFLAGS_OTHER}" )
target_include_directories(${LIBNAME} PUBLIC ${MAGMA_INCLUDE_DIRS})
target_link_directories(${LIBNAME} PUBLIC ${MAGMA_LIBRARY_DIRS})
target_link_libraries(${LIBNAME} PRIVATE ${MAGMA_LIBRARIES})
target_compile_options(${LIBNAME} PRIVATE ${MAGMA_CFLAGS_OTHER})
target_compile_options(${LIBNAME} PUBLIC -DHAVE_CUDA -DHAVE_MAGMA -DUSE_MAGMA)
#target_compile_options(${LIBNAME} PUBLIC -DHAVE_CUDA -DHAVE_MAGMA -DUSE_MAGMA)
endif(USE_MAGMA)


install(TARGETS ${LIBNAME} DESTINATION lib)
install(FILES ${INCLUDEFILES} DESTINATION include/${SRCNAME})
