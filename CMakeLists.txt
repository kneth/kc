cmake_minimum_required(VERSION 3.20)
project(kc VERSION 1.10 LANGUAGES C)

include(GNUInstallDirs)

find_package(BISON REQUIRED)
find_package(FLEX REQUIRED)

# The legacy sources still rely on pre-C99 behaviour (implicit
# declarations, tentative definitions shared between translation units).
# These flags are temporary and are removed once the code is cleaned up.
set(CMAKE_C_STANDARD 90)
set(CMAKE_C_EXTENSIONS ON)
add_compile_definitions(_PLATFORM_GCC_)
add_compile_options(-fcommon -Wno-implicit-function-declaration -Wno-implicit-int)

# Numerical library: linked into kc and into generated simulators.
add_library(kcsolvers STATIC
  Solvers/complex.c
  Solvers/eigen.c
  Solvers/integr.c
  Solvers/matrix.c
  Solvers/odeserv.c
  Solvers/odesolv.c
  Solvers/quench.c)
target_include_directories(kcsolvers PUBLIC Solvers)
target_link_libraries(kcsolvers PUBLIC m)

# main.c includes the generated parser.c, which includes the generated lex.c.
bison_target(parser src/kc.y ${CMAKE_CURRENT_BINARY_DIR}/parser.c COMPILE_FLAGS -y)
flex_target(lexer src/kc.l ${CMAKE_CURRENT_BINARY_DIR}/lex.c)
set_source_files_properties(src/main.c PROPERTIES
  OBJECT_DEPENDS "${BISON_parser_OUTPUTS};${FLEX_lexer_OUTPUTS}")

add_executable(kc
  src/codecall.c
  src/codegen.c
  src/conis.c
  src/kgode.c
  src/kncont.c
  src/main.c
  src/misc.c
  src/mixed.c
  src/symbmath.c
  src/tableman.c
  src/waves.c)
target_include_directories(kc PRIVATE src ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(kc PRIVATE kcsolvers)

configure_file(kci.in ${CMAKE_CURRENT_BINARY_DIR}/kci @ONLY)

install(TARGETS kc kcsolvers)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/kci TYPE BIN)
install(FILES
  Solvers/complex.h
  Solvers/integr.h
  Solvers/matrix.h
  Solvers/odeserv.h
  Solvers/odesolv.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/kc)
install(FILES Solvers/kksolver.c DESTINATION ${CMAKE_INSTALL_DATADIR}/kc)
