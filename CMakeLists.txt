cmake_minimum_required(VERSION 3.10)

project(yasph_proj LANGUAGES C)

set(C_STANDARD C99)
set(C_STANDARD_REQUIRED ON)

add_executable(yasph yasph.c)
target_include_directories(yasph PUBLIC "include")
target_link_libraries(yasph PRIVATE -lm -lrt -lgomp)
target_compile_options(yasph PRIVATE -Wall -Wextra -Wpedantic -O2 -fopenmp -Wno-implicit-fallthrough)

add_executable(test_io test_io.c)
target_include_directories(test_io PUBLIC "include")
target_compile_options(test_io PRIVATE -Wall -Wextra -Wpedantic)

add_executable(test_kernel test_kernel.c)
target_include_directories(test_kernel PUBLIC "include")
target_link_libraries(test_kernel PRIVATE -lm)
target_compile_options(test_kernel PRIVATE -Wall -Wextra -Wpedantic)

add_executable(test_index test_index.c)
target_include_directories(test_index PUBLIC "include")
target_link_libraries(test_index PRIVATE -lm)
target_compile_options(test_index PRIVATE -Wall -Wextra -Wpedantic -Wno-implicit-fallthrough -Wno-unused-parameter)

add_executable(test_argsio test_argsio.c)
target_include_directories(test_argsio PUBLIC "include")
target_compile_options(test_argsio PRIVATE -Wall -Wextra -Wpedantic)

add_executable(test_qtree test_qtree.c)
target_include_directories(test_qtree PUBLIC "include")
target_compile_options(test_qtree PRIVATE -Wall -Wextra -Wpedantic -Wno-unused-parameter)

enable_testing()

add_test(NAME csv-test-nonrandom COMMAND test_io csv-write-read 100 10 0 0 tempfile-nonrandom.csv)
add_test(NAME tsv-test-nonrandom COMMAND test_io tsv-write-read 80 8 0 0 tempfile-nonrandom.tsv)
add_test(NAME ssv-test-nonrandom COMMAND test_io ssv-write-read 70 7 0 0 tempfile-nonrandom.ssv)
add_test(NAME csv-test-random COMMAND test_io csv-write-read 50 5 1 0 tempfile-random.csv)
add_test(NAME tsv-test-random COMMAND test_io tsv-write-read 30 10 1 0 tempfile-random.tsv)
add_test(NAME ssv-test-random COMMAND test_io ssv-write-read 3 30 1 0 tempfile-random.ssv)
add_test(NAME csv-test-random-header COMMAND test_io csv-write-read 51 8 1 1 tempfile-random-header.csv)
add_test(NAME tsv-test-random-header COMMAND test_io tsv-write-read 241 9 1 1 tempfile-random-header.tsv)
add_test(NAME ssv-test-random-header COMMAND test_io ssv-write-read 10 3 1 1 tempfile-random-header.ssv)

set(NumAbscissas "700")
set(KernelSigmaList "0.25" "0.50" "1.00" "2.00" "4.00")
foreach(KernelSigma ${KernelSigmaList})
  add_test(NAME kernel-cubic-${KernelSigma} COMMAND test_kernel cubic ${NumAbscissas} ${KernelSigma})
  add_test(NAME kernel-quartic-${KernelSigma} COMMAND test_kernel quartic ${NumAbscissas} ${KernelSigma})
  add_test(NAME kernel-quintic-${KernelSigma} COMMAND test_kernel quintic ${NumAbscissas} ${KernelSigma})
  add_test(NAME kernel-wend2-${KernelSigma} COMMAND test_kernel wendlandc2 ${NumAbscissas} ${KernelSigma})
  add_test(NAME kernel-wend4-${KernelSigma} COMMAND test_kernel wendlandc4 ${NumAbscissas} ${KernelSigma})
  add_test(NAME kernel-wend6-${KernelSigma} COMMAND test_kernel wendlandc6 ${NumAbscissas} ${KernelSigma})
  add_test(NAME kernel-gaussian-${KernelSigma} COMMAND test_kernel gaussian ${NumAbscissas} ${KernelSigma})
  add_test(NAME kernel-sgauss-${KernelSigma} COMMAND test_kernel sgaussian ${NumAbscissas} ${KernelSigma})
endforeach()

set(NumTestPoints "5000")
add_test(NAME index-1 COMMAND test_index ${NumTestPoints} 1.0 0.0125)
add_test(NAME index-2 COMMAND test_index ${NumTestPoints} 1.0 0.025)
add_test(NAME index-3 COMMAND test_index ${NumTestPoints} 1.0 0.05)
add_test(NAME index-4 COMMAND test_index ${NumTestPoints} 1.0 0.10)
add_test(NAME index-5 COMMAND test_index ${NumTestPoints} 1.0 0.20)
add_test(NAME index-6 COMMAND test_index ${NumTestPoints} 1.0 0.40)
add_test(NAME index-7 COMMAND test_index ${NumTestPoints} 1.0 0.80)

set(QueryHalfwidth "0.05")
add_test(NAME qtree-1 COMMAND test_qtree ${NumTestPoints} 1.0 ${QueryHalfwidth} 1)  # different numbers for "max leaf"
add_test(NAME qtree-2 COMMAND test_qtree ${NumTestPoints} 1.0 ${QueryHalfwidth} 2)
add_test(NAME qtree-3 COMMAND test_qtree ${NumTestPoints} 1.0 ${QueryHalfwidth} 3)
add_test(NAME qtree-4 COMMAND test_qtree ${NumTestPoints} 1.0 ${QueryHalfwidth} 4)
add_test(NAME qtree-5 COMMAND test_qtree ${NumTestPoints} 1.0 ${QueryHalfwidth} 5)
add_test(NAME qtree-6 COMMAND test_qtree ${NumTestPoints} 1.0 ${QueryHalfwidth} 6)
add_test(NAME qtree-7 COMMAND test_qtree ${NumTestPoints} 1.0 ${QueryHalfwidth} 7)
add_test(NAME qtree-8 COMMAND test_qtree ${NumTestPoints} 1.0 ${QueryHalfwidth} 8)

set(NumTestArgsList "0" "1" "2" "3" "10" "20" "50" "100" "500" "1000")
foreach(NumTestArgs ${NumTestArgsList})
  add_test(NAME args-test-${NumTestArgs}-A COMMAND test_argsio ${NumTestArgs})
  add_test(NAME args-test-${NumTestArgs}-B COMMAND test_argsio ${NumTestArgs} --keyname)
endforeach()

add_test(NAME smoke-main-0 COMMAND yasph)
set_tests_properties(smoke-main-0 PROPERTIES WILL_FAIL true)

add_test(NAME smoke-main-1 COMMAND yasph fil1)
set_tests_properties(smoke-main-1 PROPERTIES WILL_FAIL true)

add_test(NAME smoke-main-2 COMMAND yasph fil1 fil2)
add_test(NAME smoke-main-3 COMMAND yasph fil1 fil2 a=b)

add_test(NAME box-gen-script COMMAND python3 ${PROJECT_SOURCE_DIR}/scripts/yasph-init-box.py --n 2017) 
add_test(NAME smoke-box-nosteps COMMAND yasph state.yasph params.yasph dt=1.0e-6 steps=0 threads=0) 

set(SimTestSteps "2500")
add_test(NAME sim-t1 COMMAND yasph state.yasph params.yasph dt=1.0e-6 threads=1 steps=${SimTestSteps})
add_test(NAME sim-t2 COMMAND yasph state.yasph params.yasph dt=1.0e-6 threads=2 steps=${SimTestSteps})
add_test(NAME sim-t2-logged COMMAND yasph state.yasph params.yasph dt=1.0e-6 threads=2 steps=${SimTestSteps} trace-steps=10 trace-file=simtrace-t2.yasph)
add_test(NAME log-t2-check COMMAND python3 ${PROJECT_SOURCE_DIR}/scripts/yasph-check-trace.py --trace-file simtrace-t2.yasph --trace-names energy px py cmx cmy lz --tolerance 1e-5 --force-crash --make-pdfs) 
add_test(NAME sim-t4-logged COMMAND yasph state.yasph params.yasph dt=1.0e-6 threads=4 steps=${SimTestSteps} trace-steps=10 trace-file=simtrace-t4.yasph)
add_test(NAME diff-t2-t4 COMMAND diff simtrace-t2.yasph simtrace-t4.yasph)
add_test(NAME sim-t3-snapshot COMMAND yasph state.yasph params.yasph dt=1.0e-6 threads=3 steps=${SimTestSteps} frame-steps=100 frame-file=snapshots-t3.yasph trace-steps=10 trace-file=simtrace-t3.yasph) 
add_test(NAME diff-t2-t3 COMMAND diff simtrace-t2.yasph simtrace-t3.yasph)
add_test(NAME sim-t1-snapshot COMMAND yasph state.yasph params.yasph dt=1.0e-6 threads=1 steps=${SimTestSteps} frame-steps=100 frame-file=snapshots-t1.yasph trace-steps=10 trace-file=simtrace-t1.yasph) 
add_test(NAME diff-t2-t1 COMMAND diff simtrace-t2.yasph simtrace-t1.yasph) 
add_test(NAME cmp-t1-t3 COMMAND cmp snapshots-t1.yasph snapshots-t3.yasph)  

add_test(NAME barrier-gen-script COMMAND python3 ${PROJECT_SOURCE_DIR}/scripts/yasph-init-box.py --n 1999 --particle-file state2.yasph --parameter-file params2.yasph --barrier-radius 2.0) 
add_test(NAME smoke-barrier-nosteps COMMAND yasph state2.yasph params2.yasph dt=1.0e-6 steps=0) 

add_test(NAME sim-t2-logged-quad COMMAND yasph state.yasph params.yasph index-name=quad dt=1.0e-6 threads=2 steps=${SimTestSteps} trace-steps=10 trace-file=simtrace-t2-quad.yasph)
add_test(NAME quad-check-1 COMMAND python3 ${PROJECT_SOURCE_DIR}/scripts/yasph-cmp-traces.py --trace-1 simtrace-t2.yasph --trace-2 simtrace-t2-quad.yasph)
