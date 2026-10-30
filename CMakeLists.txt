cmake_minimum_required(VERSION 3.20)
project(kc VERSION 1.10 LANGUAGES C)

include(GNUInstallDirs)

find_package(BISON REQUIRED)
find_package(FLEX REQUIRED)

if(NOT CMAKE_C_STANDARD)
  set(CMAKE_C_STANDARD 17)
endif()
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)
add_compile_definitions(_POSIX_C_SOURCE=200809L)  # getopt, fileno
add_compile_options(-Wall -Wextra -Wpedantic)
# No fused multiply-add: mode 1 results must not depend on compiler or CPU.
add_compile_options(-ffp-contract=off)

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

# Golden-output tests: one test per example model and code generator mode.
enable_testing()
add_executable(runner test/runner.c)
set_target_properties(runner PROPERTIES C_STANDARD 17 C_EXTENSIONS OFF)
target_link_libraries(runner PRIVATE m)
target_compile_definitions(runner PRIVATE
  KC_CC="${CMAKE_C_COMPILER}"
  KC_SIM_CFLAGS="-w -std=gnu89"
  KC_SOLVERS="${CMAKE_CURRENT_SOURCE_DIR}/Solvers"
  KC_SOLVERS_LIB="$<TARGET_FILE:kcsolvers>"
  KC_SKIP_LIST="${CMAKE_CURRENT_SOURCE_DIR}/test/skip.txt")
add_dependencies(runner kcsolvers)

file(GLOB models LIST_DIRECTORIES false ${CMAKE_CURRENT_SOURCE_DIR}/test/*)
list(FILTER models EXCLUDE REGEX "\\.(c|txt)$")
set(test_work ${CMAKE_CURRENT_BINARY_DIR}/test-work)
file(MAKE_DIRECTORY ${test_work})
foreach(model IN LISTS models)
  get_filename_component(name ${model} NAME)
  foreach(mode 1 2 3 4 5)
    add_test(NAME ${name}.m${mode}
      COMMAND runner $<TARGET_FILE:kc> ${model} ${mode} ${CMAKE_CURRENT_SOURCE_DIR}/test/expected
      WORKING_DIRECTORY ${test_work})
    set_tests_properties(${name}.m${mode} PROPERTIES SKIP_RETURN_CODE 77 TIMEOUT 120)
  endforeach()
endforeach()
