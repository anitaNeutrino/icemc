cmake_minimum_required(VERSION 2.8.10)
project(icemc)

set(libname "icemclib")

### Site-specific options

SET(BOOST_INCLUDE_DIR "" CACHE String "BOOST include path (if different from system path)")  
#option(USE_CPP11  "Use the C++11 standard" ON)
option(ANITA3_EVENTREADER "Load ANITA3 EventReader" ON) 


## ROOT 
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS MathMore MathCore RIO Hist Tree Net)

message("ROOT_VERSION is set to ${ROOT_VERSION}")
include(${ROOT_USE_FILE})

include_directories(${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS} ${FFTW_INCLUDES} ${ZLIB_INCLUDES} ${UTIL_INC_DIR})
add_definitions(${ROOT_CXX_FLAGS})
#include_directories(${PROJECT_SOURCE_DIR})

message("cmake_cxx_flags = ${CMAKE_CXX_FLAGS}")

string(FIND ${CMAKE_CXX_FLAGS} "-std=" NEED_SET_CPP_FLAGS)
if(${NEED_SET_CPP_FLAGS} EQUAL -1)
  if(CMAKE_COMPILER_IS_GNUCC)
    message(${CMAKE_CXX_COMPILER_VERSION})
    if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 4.5.0 OR CMAKE_CXX_COMPILER_VERSION VERSION_EQUAL 4.5.0)
      message(" -std=c++11 ") 
      add_definitions(-std=c++11)
    else()
      ##DRAGONS BE HERE
      add_definitions(-Dconstexpr=const)
      message(" -std=c++0x ")
      add_definitions(-std=c++0x)
    endif()
  endif()
endif()


#string(FIND ${CMAKE_CXX_FLAGS} "-std=c++1" USE_MODERN_CPP)
#if(${USE_MODERN_CPP} EQUAL -1)
#  ##DRAGONS BE HERE
#  add_definitions(-Dconstexpr=const)
#endif()


set(DICTNAME G__${libname})

#set(DICT classdict)
set(CLASS_HEADERS rx.hpp Taumodel.hh Settings.h blazars/fava.h icemc_random.h)

ROOT_GENERATE_DICTIONARY(${DICTNAME} ${CLASS_HEADERS} LINKDEF LinkDef.h)
#ROOT_GENERATE_DICTIONARY(${DICT} ${CLASS_HEADERS} LINKDEF LinkDef.h)
FILE (GLOB SOURCE_FILES "*.cc") 

SET(SOURCE vector.cc position.cc earthmodel.cc balloon.cc icemodel.cc signal.cc
  ray.cc Spectra.cc anita.cc roughness.cc secondaries.cc Primaries.cc Tools.cc
  counting.cc Settings.cc Taumodel.cc screen.cc GlobalTrigger.cc
  ChanTrigger.cc SimulatedSignal.cc EnvironmentVariable.cc source.cc random.cc)

#SET(SOURCES vector.cc position.cc earthmodel.cc balloon.cc icemodel.cc signal.cc
#  ray.cc Spectra.cc anita.cc roughness.cc secondaries.cc Primaries.cc Tools.cc
#  counting.cc Settings.cc Taumodel.cc screen.cc GlobalTrigger.cc
#  ChanTrigger.cc SimulatedSignal.cc EnvironmentVariable.cc source.cc random.cc)


if( ${ROOT_VERSION} VERSION_GREATER "5.99.99")
    add_custom_target(${DICTNAME}.pcm DEPENDS ${DICT})
#    add_custom_target(${DICT}.pcm DEPENDS ${DICT})
endif()


set(LIBS ${ROOT_LIBRARIES} MathMore ) 

if (ANITA3_EVENTREADER) 
  add_definitions(-DANITA3_EVENTREADER)
  set(LIBS ${LIBS} AnitaEvent )
  add_definitions(-DANITA3_EVENTCORRELATOR)
  set(LIBS ${LIBS} AnitaCorrelator )
endif()


if(DEFINED ENV{ANITA_UTIL_INSTALL_DIR})
    message("ANITA_UTIL_INSTALL_DIR is found and set to $ENV{ANITA_UTIL_INSTALL_DIR}")
    add_definitions(-DANITA_UTIL_EXISTS)
    set (LIBS ${LIBS} RootFftwWrapper ) 
    set(ANITA_UTIL_LIB_DIR $ENV{ANITA_UTIL_INSTALL_DIR}/lib)
    set(ANITA_UTIL_INC_DIR $ENV{ANITA_UTIL_INSTALL_DIR}/include)
    set(LD_ANITA_UTIL $ENV{ANITA_UTIL_INSTALL_DIR}/lib)
    set(INC_ANITA_UTIL $ENV{ANITA_UTIL_INSTALL_DIR}/include)
    set(ANITA_UTIL_ETC_DIR $ENV{ANITA_UTIL_INSTALL_DIR}/etc)
    link_directories(${ANITA_UTIL_LIB_DIR})
    include_directories(${INC_ANITA_UTIL})
 endif()


set(CMAKE_CXX_FLAGS_DEFAULT 
  "-O2 -g -pipe -m64 -pthread -W -Wall -Wextra -Woverloaded-virtual" 
  CACHE STRING "c++ Flags used during default icemc builds" 
  FORCE ) 

mark_as_advanced ( CMAKE_CXX_FLAGS_DEFAULT CMAKE_C_FLAGS_DEFAULT CMAKE_EXE_LINKER_FLAGS_DEFAULT CMAKE_SHARED_LINKER_FLAGS_DEFAULT) 

if (NOT CMAKE_BUILD_TYPE) 
  set (CMAKE_BUILD_TYPE Default
    CACHE STRING "Choose type of build: None Debug Release RelWithDebInfo MinSizeRel Default"
    FORCE ) 
endif()

option (USE_HEALPIX "Enable Healpix Segmentation, requires healpix to be
installed on your system in a place that can be found easily." OFF)

option (VECTORIZE "Enable manual vectorization, does nothing if ANITA_UTIL is
not available." OFF)


if (USE_HEALPIX)
  find_package(PkgConfig) 
  pkg_search_module(HEALPIX REQUIRED healpix_cxx) 
  add_definitions  (-DUSE_HEALPIX )
  SET(LIBS ${LIBS} ${HEALPIX_LIBRARIES})
  include_directories(${HEALPIX_INCLUDE_DIRS})
  link_directories(${HEALPIX_LIBRARY_DIRS})
endif()

if (VECTORIZE) 
  add_definitions(-DVECTORIZE) 
  if(CMAKE_COMPILER_IS_GNUCXX) 
      ### someone should do this for clang if they want it to be as fast as possible 
      if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 4.1)
        if(NATIVE_ARCH) 
          add_definitions( -march=native)
        endif()
      endif()
      add_definitions(  -fabi-version=0 ) 
      include_directories(${INC_ANITA_UTIL}/vectorclass) 
  endif() 
endif() 

#---Create a shared library with generated dictionary
add_library(${libname} SHARED ${SOURCE} ${DICTNAME}.cxx)
target_link_libraries(${libname} ${ZLIB_LIBRARIES} RootFftwWrapper AnitaEvent AnitaCorrelator AnitaAnalysis ${ROOT_LIBRARIES} ${FFTW_LIBRARIES})

##We need to compile stuff into a library unless we want it to be compiled a
##million times over 
#add_library(icemc_objs STATIC ${SOURCES} ${DICT}.cxx)

macro(add_binary bin) 
  add_executable(${bin} ${bin}.cc)
  target_link_libraries(${bin} ${libname} ${LIBS})
#  target_link_libraries(${bin} icemc_objs ${LIBS})
endmacro() 

add_binary(icemc) 
add_binary(testTrigger) 
add_binary(testSettings) 
add_binary(testEAS) 
add_binary(testWAIS) 
add_binary(testThermalNoise)
add_binary(testInputAfterAntenna)
