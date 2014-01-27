# TODO add partitiontree target
cmake_minimum_required( VERSION 2.8 )

project( mop )

IF(NOT CMAKE_BUILD_TYPE)
	SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: Debug Release Profile" FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wextra -Wno-long-long -Wno-unknown-pragmas -g -ggdb -O0 -DDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-Wall -Wextra -Wno-long-long -Wno-unknown-pragmas -O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_PROFILE "-Wall -Wextra -Wno-long-long -Wno-unknown-pragmas -O3 -g -ggdb -DNDEBUG")

set( Boost_USE_STATIC_LIBS ON )
find_package( Boost 1.46.0 REQUIRED )

set( common_src
	src/verbose.cpp
)

set( common_hdr
	src/molecule.h
	src/verbose.h
)

set ( bk_hdr
	src/bronkerbosch.h
)

set ( bk_src
	src/bk.cpp
)

set( cgp_hdr
	src/cgp/solver.h
	src/cgp/solvertreewidthdp.h
	src/treedecomposition.h
)

set( cgp_src
	src/partition.cpp
)

set( sym_hdr
	src/symmetrization/symmetrization.h
	src/symmetrization/symmetrizationautomorphism.h
)

set( sym_src
	src/symmetrize.cpp
  #	../GNA/webserver.cpp
)

set( atb2lgf_src
	src/atb2lgf.cpp
)

set( atb2lgf_hdr
	src/cgp/solverdummy.h
)

set( readsolution_src
	src/readsolution.cpp
)

set( readsolution_hdr
	src/cgp/solverdummy.h
)

set( fragments_src
	src/fragments.cpp
)

set( fragments_hdr
	src/fragments/product.h
	src/bronkerbosch.h
	src/fragments/bronkerboschconnected.h
)

set( LIBLEMON_ROOT "$ENV{HOME}/lemon" CACHE PATH "Additional search directory for lemon library" )

set( CommonLibs
        m
        emon
)

if( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libstdc++ " )
endif()

configure_file( src/config.h.in ${PROJECT_SOURCE_DIR}/src/config.h )
include_directories( "src" "${LIBLEMON_ROOT}/include" )
link_directories( "${LIBLEMON_ROOT}/lib" )

add_executable( partition ${common_src} ${common_hdr} ${cgp_src} ${cgp_hdr} ) 
target_link_libraries( partition ${CommonLibs} )

#add_executable( symmetrize ${common_src} ${common_hdr} ${sym_src} ${sym_hdr} )

add_executable( readsolution ${common_src} ${common_hdr} ${readsolution_src} ${readsolution_hdr} )
target_link_libraries( readsolution ${CommonLibs} )

add_executable( atb2lgf ${common_src} ${common_hdr} ${atb2lgf_src} ${atb2lgf_hdr} )
target_link_libraries( atb2lgf ${CommonLibs} )

add_executable( fragments ${common_src} ${common_hdr} ${fragments_src} ${fragments_hdr} )
target_link_libraries( fragments ${CommonLibs} )

add_executable( bk ${common_src} ${common_hdr} ${bk_src} ${bk_hdr} )
target_link_libraries( bk ${CommonLibs} )
