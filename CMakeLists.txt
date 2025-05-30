# For Debian currently with
#
#   cd build
#   cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo -DZIG=OFF ..
#   cmake --build .
#   ctest .
#   cmake --install .
#
cmake_minimum_required(VERSION 3.16)
project(vcflib)

enable_testing()

set(CMAKE_CXX_STANDARD 17)
include(ExternalProject)
include(FeatureSummary)
include(GNUInstallDirs)

find_package(BZip2 REQUIRED)   # htslib
find_package(LibLZMA REQUIRED) # htslib
find_package(PkgConfig REQUIRED)
find_package(ZLIB REQUIRED)
find_package(pybind11 CONFIG)

include(GNUInstallDirs)

feature_summary(
  FATAL_ON_MISSING_REQUIRED_PACKAGES
  WHAT REQUIRED_PACKAGES_NOT_FOUND)

if (BUILD_STATIC)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  set(BUILD_SHARED_LIBS OFF)
  set(CMAKE_EXE_LINKER_FLAGS "-static")
  set(CMAKE_POSITION_INDEPENDENT_CODE OFF)
else ()
  set(CMAKE_POSITION_INDEPENDENT_CODE ON) # for pybind11
endif ()

# ---- Options

option(BUILD_STATIC "Build static binary" OFF)
option(BUILD_DOC "Build documentation" ON)
option(OPENMP "Enable OpenMP" ON) # disabling does not work because of vcfwave
option(PROFILING "Enable profiling" OFF)
option(GPROF "Enable gprof profiling" OFF)
option(ASAN "Use address sanitiser" OFF)
option(ZIG "Set to OFF to disable the zig code" ON)
option(WFA_GITMODULE "Force local git submodule for WFA2LIB" OFF) # disabled by default now
include(CheckIPOSupported) # adds lto
check_ipo_supported(RESULT ipo_supported OUTPUT output)

set(EXTRA_FLAGS " ") # otherwise it injects OFF

# ---- Dependencies

if(OPENMP)
  include(FindOpenMP)
endif(OPENMP)

find_package(ZLIB)
set_package_properties(ZLIB PROPERTIES TYPE REQUIRED)
find_package(Threads)
set_package_properties(Threads PROPERTIES TYPE REQUIRED)

pkg_check_modules(htslib IMPORTED_TARGET htslib)   # Optionally builds from contrib/
pkg_check_modules(tabixpp IMPORTED_TARGET tabixpp) # Optionally builds from contrib/

# ---- Build switches
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ${ipo_supported})

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING
          "Choose the type of build, options are: Release|Debug|RelWithDebInfo (for distros)." FORCE)
endif()

if (${CMAKE_BUILD_TYPE} MATCHES Release)
  set(EXTRA_FLAGS "-march=native -D_FILE_OFFSET_BITS=64")
  # set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG") # reset CXX_FLAGS to replace -O3 with -Ofast
endif()

if ((${CMAKE_BUILD_TYPE} MATCHES Release) OR (${CMAKE_BUILD_TYPE} MATCHES RelWithDebInfo))
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${EXTRA_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${EXTRA_FLAGS}")
endif ()

if (${CMAKE_BUILD_TYPE} MATCHES "Debug")
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${EXTRA_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${EXTRA_FLAGS}")
  add_definitions(-Wfatal-errors) # stop after first error
endif ()

if (ASAN)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=address  -fno-omit-frame-pointer -fno-common")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address  -fno-omit-frame-pointer -fno-common")
endif(ASAN)

if(PROFILING)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
endif(PROFILING)

if(GPROF)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
endif(GPROF)

if (ZIG)
  find_program(ZIG_EXE NAMES "zig")
  if (NOT ZIG_EXE)
    MESSAGE(FATAL_ERROR "zig binary not found in PATH. zig is used for vcfcreatemulti's latest features. Either use cmake -DZIG=OFF option or add zig to the PATH")
  endif (NOT ZIG_EXE)
endif(ZIG)

# ---- Include files

include_directories(include)
include_directories(contrib/fastahack)
include_directories(contrib/intervaltree)
include_directories(contrib/smithwaterman)
include_directories(contrib/filevercmp)
include_directories(contrib/c-progress-bar)

if(NOT htslib_FOUND)
  message(STATUS "Using included htslib")
  include(FindCURL) # for htslib
  set(htslib_LOCAL contrib/tabixpp/htslib)
  set(tabixpp_FOUND OFF) # also build tabixpp if htslib is missing
endif()

if(NOT tabixpp_FOUND)
  message(STATUS "Using included tabixpp")
  set(tabixpp_LOCAL contrib/tabixpp)
  include_directories(contrib/tabixpp)
  set(tabixpp_SOURCE
    contrib/tabixpp/tabix.cpp
  )
endif()

if(WFA_GITMODULE)
  set(WFA_LOCAL contrib/WFA2-lib)
endif()

file(GLOB INCLUDES
  src/*.h*
  contrib/intervaltree/*.h*
  contrib/smithwaterman/*.h*
  contrib/fastahack/*.h*
  contrib/filevercmp/*.h*
  )

set(vcfwfa_SOURCE
    src/legacy.cpp # introduces a WFA dependency
    src/vcf-wfa.cpp
)

set(vcflib_SOURCE
    ${vcfwfa_SOURCE}
    src/vcf-c-api.cpp
    src/Variant.cpp
    # src/canonicalize.cpp --- merged back into Variant.h, see https://github.com/vcflib/vcflib/issues/417
    src/rnglib.cpp
    src/var.cpp
    src/pdflib.cpp
    src/cdflib.cpp
    src/split.cpp
    src/stats.cpp
    src/index.cpp
    src/phase.cpp
    src/rkmh.cpp
    src/murmur3.cpp
    src/LeftAlign.cpp
    src/cigar.cpp
    src/allele.cpp
    contrib/fastahack/Fasta.cpp
    contrib/smithwaterman/SmithWatermanGotoh.cpp
    contrib/smithwaterman/Repeats.cpp
    contrib/smithwaterman/IndelAllele.cpp
    contrib/smithwaterman/disorder.cpp
    contrib/smithwaterman/LeftAlign.cpp
    contrib/fsom/fsom.c
    contrib/filevercmp/filevercmp.c
    contrib/c-progress-bar/progress.c
)

if (tabixpp_LOCAL) # add the tabixpp source file
    list(APPEND vcflib_SOURCE ${tabixpp_SOURCE})
    list(APPEND INCLUDES ${TABIXPP_LOCAL}/tabix.hpp)
endif()

add_library(vcflib
    ${vcflib_SOURCE}
    )

set(BINS
    vcfecho
    dumpContigsFromHeader
    bFst
    pVst
    hapLrt
    popStats
    wcFst
    iHS
    segmentFst
    segmentIhs
    genotypeSummary
    sequenceDiversity
    pFst
    smoother
    vcfld
    plotHaps
    abba-baba
    permuteGPAT++
    permuteSmooth
    normalize-iHS
    meltEHH
    vcfnullgenofields
    vcfaltcount
    vcfhetcount
    vcfhethomratio
    vcffilter
    vcf2tsv
    vcfgenotypes
    vcfannotategenotypes
    vcfcommonsamples
    vcfremovesamples
    vcfkeepsamples
    vcfsamplenames
    vcfgenotypecompare
    vcffixup
    vcfclassify
    vcfsamplediff
    vcfremoveaberrantgenotypes
    vcfrandom
    vcfflatten
    vcfprimers
    vcfnumalt
    vcfintersect
    vcfannotate
    vcfoverlay
    vcfaddinfo
    vcfkeepinfo
    vcfkeepgeno
    vcfafpath
    vcfcountalleles
    vcflength
    vcfdistance
    vcfrandomsample
    vcfentropy
    vcfglxgt
    vcfcheck
    vcfstreamsort
    vcfuniq
    vcfuniqalleles
    vcfremap
    vcf2fasta
    vcfsitesummarize
    vcfbreakmulti
    vcfevenregions
    vcfcat
    vcfgenosummarize
    vcfgenosamplenames
    vcfgeno2haplo
    vcfleftalign
    vcfcombine
    vcfgeno2alleles
    vcfindex
    vcf2dag
    vcfsample2info
    vcfqual2info
    vcfinfo2qual
    vcfglbound
    vcfinfosummarize
    vcfcreatemulti
)

set(WFBINS # Introduce wavefront dependency
    vcfallelicprimitives
    vcfcleancomplex
    vcfparsealts
    vcfroc
    vcfstats
    vcfwave
  )

set(SCRIPTS
    bed2region
    bgziptabix
    vcf2bed.py
    vcf2sqlite.py
    vcfbiallelic
    vcfclearid
    vcfclearinfo
    vcfcomplex
    vcffirstheader
    vcfgtcompare.sh
    vcfindelproximity
    vcfindels
    vcfjoincalls
    vcfmultiallelic
    vcfmultiway
    vcfmultiwayscripts
    vcfnobiallelicsnps
    vcfnoindels
    vcfnosnps
    vcfnulldotslashdot
    vcfplotaltdiscrepancy.r
    vcfplotaltdiscrepancy.sh
    vcfplotsitediscrepancy.r
    vcfplottstv.sh
    vcfprintaltdiscrepancy.r
    vcfprintaltdiscrepancy.sh
    vcfqualfilter
    vcfregionreduce
    vcfregionreduce_and_cut
    vcfregionreduce_pipe
    vcfregionreduce_uncompressed
    vcfremovenonATGC
    vcfsnps
    vcfsort
    vcf_strip_extra_headers
    vcfvarstats
    )

# ---- Get version

file (STRINGS "VERSION" BUILD_NUMBER)
add_definitions(-DVCFLIB_VERSION="${BUILD_NUMBER}")
add_definitions(-DVERSION="${BUILD_NUMBER}")

# ---- Build htslib
#
# Note by default we use the distributed htslib! These are
# the older instructions which allow for a local build from
# git submodules

if (htslib_LOCAL)

  include_directories(${htslib_LOCAL})

  set(flags "-O2 -g -fPIC")
  ExternalProject_Add(htslib-EXT
    SOURCE_DIR "${CMAKE_SOURCE_DIR}/${htslib_LOCAL}"
    UPDATE_COMMAND autoreconf -i
    CONFIGURE_COMMAND ./configure --disable-s3
    INSTALL_COMMAND ""
    BUILD_IN_SOURCE ON
    BUILD_COMMAND $(MAKE) CFLAGS=${flags} lib-static
  )
  ExternalProject_Get_property(htslib-EXT SOURCE_DIR)
  ExternalProject_Get_property(htslib-EXT INSTALL_DIR)
  set(htslib_INCLUDE "${SOURCE_DIR}")
  set(htslib_static "${SOURCE_DIR}/libhts.a")
  set(htslib_lib "${htslib_static}")

  add_library(htslib STATIC IMPORTED)
  set_target_properties(htslib PROPERTIES IMPORTED_LOCATION ${htslib_lib})
  set_target_properties(htslib PROPERTIES INCLUDE_DIRECTORIES ${htslib_INCLUDE})
  add_dependencies(htslib htslib-EXT)

# If the user wants to configure our HTSlib to build with libddeflate, we need
# to make sure to link against libdeflate as a transitive dependency. To do
# that, pass -Dhtslib_EXTRA_LIBS="-ldeflate" when configuring the project with
# cmake.
# TODO: Stop vendoring in htslib and just use find_package

  # set(htslib_EXTRA_LIBS "-lcurl" CACHE STRING "Library flags needed to link with htslib's dependencies, for chosen configuration")
  # set_property(TARGET htslib PROPERTY INTERFACE_LINK_LIBRARIES ${htslib_EXTRA_LIBS})

endif(htslib_LOCAL)

if(WFA_GITMODULE)
  message(STATUS "Using included libwfa")
  set(WFA_INCLUDE_DIRS ${WFA_LOCAL})
  add_subdirectory(${WFA_LOCAL})
  set(WFALIB wfa2) # pick up the wfa2 lib target from the included CMakeLists.txt
else(WFA_GITMODULE)
  find_path(WFA_INCLUDE_DIRS wavefront/wfa.hpp PATH_SUFFIXES wfa2lib)
  find_library(WFALIB wfa2 wfa REQUIRED) # distro search for shared lib
endif(WFA_GITMODULE)

# if(NOT WFALIB)
#   message(STATUS "ERROR: Can not find wfa2lib! Make sure it is installed or use the git submodule instead")
# endif()

include_directories(${WFA_INCLUDE_DIRS})
MESSAGE(STATUS "WFA using include ${WFA_INCLUDE_DIRS}")

# ZIG VCF imported library is part of VCFLIB source tree. We designate it an
# external project so we can run 'zig build'

if (ZIG)
ExternalProject_Add(ZIG-EXT
    SOURCE_DIR "${CMAKE_SOURCE_DIR}/src/zig"
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    INSTALL_COMMAND ""
    BUILD_IN_SOURCE ON
    BUILD_ALWAYS ON
    BUILD_COMMAND ${CMAKE_SOURCE_DIR}/src/zig/compile.sh # -Drelease-fast=true -freference-trace
)
ExternalProject_Get_property(ZIG-EXT SOURCE_DIR)
set(ZIG_INCLUDE_DIRS ${SOURCE_DIR})
set(ZIG_LINK_LIBRARIES ${SOURCE_DIR}/zig-out/lib/libzig.a)

add_library(ZIGCPP_STATIC_LIB STATIC IMPORTED)
set_target_properties(ZIGCPP_STATIC_LIB PROPERTIES IMPORTED_LOCATION ${ZIG_LINK_LIBRARIES})
add_dependencies(ZIGCPP_STATIC_LIB ZIG-EXT vcf.zig)

include_directories(${ZIG_INCLUDE_DIRS})

else (ZIG)
  set (ZIG_LINK_LIBRARIES )
  add_definitions(-DNO_ZIG=1)
endif (ZIG)

if (htslib_LOCAL)
  target_link_libraries(vcflib htslib)
  set(vcflib_LIBS curl deflate) # local build requires linking curl explicitly
else ()
  target_link_libraries(vcflib PkgConfig::htslib)
  set(vcflib_LIBS hts)
endif()

list(APPEND htslib_dependencies
  lzma
  z
  bz2
)

list(APPEND vcflib_LIBS
  ${CMAKE_THREAD_LIBS_INIT}
  ${htslib_dependencies} # linking htslib for everything is not needed if we have our own bgzip reader
)

#  # ${WFA_LINK_LIBRARIES}


if (NOT tabixpp_LOCAL)
  target_link_libraries(vcflib PkgConfig::tabixpp)
endif()

if(OPENMP)
  target_link_libraries(vcflib OpenMP::OpenMP_CXX)
endif()

# ---- Build all

if (BUILD_SHARED_LIBS)
  target_link_libraries(vcflib PUBLIC ${vcflib_LIBS} ${WFALIB})
endif()

if (NOT BUILD_ONLY_LIB)
  foreach(BIN ${BINS})
    add_executable(${BIN} src/${BIN}.cpp)
    add_dependencies(${BIN} vcflib)
    target_link_libraries(${BIN} PUBLIC ${vcflib_LIBS} vcflib)
  endforeach(BIN ${BINS})
  foreach(WFBIN ${WFBINS})
    add_executable(${WFBIN} src/${WFBIN}.cpp ${vcfwfa_SOURCE})
    add_dependencies(${WFBIN} vcflib)
    target_link_libraries(${WFBIN} PUBLIC ${vcflib_LIBS} vcflib)
    target_link_libraries(${WFBIN} PUBLIC ${vcflib_LIBS} ${WFALIB})
  endforeach(WFBIN ${BINS})
  if (ZIG)
    # list(APPEND vcflib_LIBS ${ZIG_LINK_LIBRARIES})
    target_link_libraries(vcfcreatemulti PUBLIC ${vcflib_LIBS} ${ZIG_LINK_LIBRARIES} vcflib)
  endif()
  install(TARGETS ${BINS} ${WFBINS} RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

  # ---- Copy scripts
  foreach(SCRIPT ${SCRIPTS})
    install(PROGRAMS ./scripts/${SCRIPT} DESTINATION ${CMAKE_INSTALL_BINDIR} RENAME ${SCRIPT})
  endforeach(SCRIPT ${SCRIPTS})
endif()


if (NOT BUILD_STATIC) # only when not building static
  # ---- Python bindings - mostly for testing at this stage
  pybind11_add_module(pyvcflib "${CMAKE_SOURCE_DIR}/src/pythonffi.cpp")
  # add_dependencies(pyvcflib CURL::libcurl)
  target_link_libraries(pyvcflib PUBLIC vcflib ${vcflib_LIBS} ${WFALIB})
  install(TARGETS pyvcflib LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

function(add_pytest TEST_FILE)
  add_test(
      NAME ${TEST_FILE}
      COMMAND python3 tests/${TEST_FILE}.py
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test
  )
  SET_TESTS_PROPERTIES(${TEST_FILE}
      PROPERTIES ENVIRONMENT "LSAN_OPTIONS=suppressions=../test/libasan-suppress.txt;PYTHONPATH=${PROJECT_SOURCE_DIR}/build;LD_LIBRARY_PATH=$ENV{LIBRARY_PATH}:$ENV{GUIX_PROFILE}/lib;PATH=${CMAKE_BINARY_DIR}/bin:$ENV{PATH}")
endfunction()

function(add_pydoctest TEST_FILE)
  add_test(
      NAME ${TEST_FILE}
      COMMAND python3 -m doctest -o NORMALIZE_WHITESPACE -o REPORT_UDIFF pytest/${TEST_FILE}.md
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test
  )
  SET_TESTS_PROPERTIES(${TEST_FILE}
      # PROPERTIES ENVIRONMENT "ASAN_OPTIONS=detect_leaks=1:symbolize=1;LSAN_OPTIONS=verbosity=2:log_threads=1:suppressions=../test/libasan-suppress.txt")
      PROPERTIES ENVIRONMENT "LSAN_OPTIONS=suppressions=../test/libasan-suppress.txt;PYTHONPATH=${PROJECT_SOURCE_DIR}/build:${PROJECT_SOURCE_DIR}/test/pytest;LD_LIBRARY_PATH=$ENV{LIBRARY_PATH}:$ENV{GUIX_PROFILE}/lib;PATH=${CMAKE_BINARY_DIR}:$ENV{PATH}")
endfunction()

function(add_doctest TEST_FILE)
  add_test(
      NAME ${TEST_FILE}
      COMMAND python3 -m doctest -o NORMALIZE_WHITESPACE -o REPORT_UDIFF ../${TEST_FILE}.md
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test
  )
endfunction()

function(add_pydoctest_fullname TEST_FILE)
  add_test(
      NAME ${TEST_FILE}
      COMMAND python3 -m doctest -o NORMALIZE_WHITESPACE -o REPORT_UDIFF ${TEST_FILE}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test
  )
endfunction()

if (ZIG_TEST) # currently disabled
function(add_zigtest)
  add_test(
    NAME zigvcf
    COMMAND zig build test
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/src/zig
    )
endfunction()

  add_zigtest()
endif(ZIG_TEST)

add_pytest(realign)
add_pydoctest(pyvcflib)
add_pydoctest(vcflib-api)
add_pydoctest(vcf2tsv)
add_pydoctest(vcfallelicprimitives)
add_pydoctest(vcfwave)
add_pydoctest(vcffilter)
if (ZIG)
  add_pydoctest(vcfcreatemulti)
endif (ZIG)

add_pydoctest(vcfnulldotslashdot)
add_doctest(doc/vcfintersect)

# ---- Build documentation
#
# Generates man pages for the python doctests. Don't need
# to run every time so it is a separate command. For pandoc logic see
# https://www.howtogeek.com/682871/how-to-create-a-man-page-on-linux/
#
# cmake --build . --target man ; cmake --install .

find_program(PANDOC pandoc)

set(pytest vcfnulldotslashdot vcfallelicprimitives vcfwave)

# vcfcreatemulti)

if (PANDOC)
    # note the option ALL which allows to build the docs together with the application
    add_custom_target( man
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMAND ruby ./test/scripts/bin2md.rb --man test/scripts/bin2md-template.erb
        COMMAND ruby ./test/scripts/bin2md.rb --man --index
        # overwrite markdown from tests
        COMMAND ${CMAKE_COMMAND} -E copy ./test/pytest/vcfnulldotslashdot.md ./doc/vcfnulldotslashdot.md
        COMMAND ${CMAKE_COMMAND} -E copy ./test/pytest/vcfallelicprimitives.md ./doc/vcfallelicprimitives.md
        COMMAND ${CMAKE_COMMAND} -E copy ./test/pytest/vcfwave.md ./doc/vcfwave.md
        COMMAND ${CMAKE_COMMAND} -E copy ./test/pytest/vcfcreatemulti.md ./doc/vcfcreatemulti.md
        COMMAND ${CMAKE_COMMAND} -E copy ./test/pytest/vcflib-api.md ./doc/vcflib-api.md
        COMMAND ${CMAKE_COMMAND} -E copy ./test/pytest/pyvcflib.md ./doc/pyvcflib.md
        # generate man pages from markdown
        COMMAND ruby ./test/scripts/md2man
        # regenerate to allow for URLs in markdown docs
        COMMAND ruby ./test/scripts/bin2md.rb test/scripts/bin2md-template.erb
        COMMAND ruby ./test/scripts/bin2md.rb --index
        COMMAND ${CMAKE_COMMAND} -E copy ./test/pytest/vcfnulldotslashdot.md ./doc/vcfnulldotslashdot.md
        COMMAND ${CMAKE_COMMAND} -E copy ./test/pytest/vcfallelicprimitives.md ./doc/vcfallelicprimitives.md
        COMMAND ${CMAKE_COMMAND} -E copy ./test/pytest/vcfwave.md ./doc/vcfwave.md
        COMMAND ${CMAKE_COMMAND} -E copy ./test/pytest/vcfcreatemulti.md ./doc/vcfcreatemulti.md
        # Documented test cases for Python bindings
        COMMAND ${CMAKE_COMMAND} -E copy ./test/pytest/pyvcflib.md ./doc/pyvcflib.md
        COMMAND ${CMAKE_COMMAND} -E copy ./test/pytest/vcflib-api.md ./doc/vcflib-api.md
    )
else (PANDOC)
  message("Pandoc needs to be installed to generate the man pages")
endif (PANDOC)

# ---- Install

install(TARGETS vcflib ARCHIVE DESTINATION ${CMAKE_INSTALL_BINDIR})
if(WFA_GITMODULE)
  install(TARGETS ${WFALIB} ARCHIVE DESTINATION ${CMAKE_INSTALL_BINDIR})
endif(WFA_GITMODULE)

install(FILES ${INCLUDES} DESTINATION include/vcflib)

install(DIRECTORY ${CMAKE_SOURCE_DIR}/man/ DESTINATION ${CMAKE_INSTALL_MANDIR}/man1)
