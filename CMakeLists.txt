cmake_minimum_required(VERSION 2.8.10)
project(icemc)

### Site-specific options

SET(BOOST_INCLUDE_DIR "" CACHE String "BOOST include path (if different from system path)")  
#option(USE_CPP11  "Use the C++11 standard" ON)
option(ANITA3_EVENTREADER "Load ANITA3 EventReader" ON) 


## ROOT 
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS MathMore MathCore RIO Hist Tree Net)

message("ROOT_VERSION is set to ${ROOT_VERSION}")
include(${ROOT_USE_FILE})

include_directories(${PROJECT_SOURCE_DIR})

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




set(DICT classdict)
set(CLASS_HEADERS rx.hpp Taumodel.hh)

ROOT_GENERATE_DICTIONARY(${DICT} ${CLASS_HEADERS} LINKDEF LinkDef.h)
FILE (GLOB SOURCE_FILES "*.cc") 

SET(SOURCES vector.cc position.cc earthmodel.cc balloon.cc icemodel.cc signal.cc
  ray.cc Spectra.cc anita.cc roughness.cc secondaries.cc Primaries.cc Tools.cc
  counting.cc Settings.cc Taumodel.cc screen.cc GlobalTrigger.cc
  ChanTrigger.cc SimulatedSignal.cc)


if( ${ROOT_VERSION} VERSION_GREATER "5.99.99")
    add_custom_target(${DICT}.pcm DEPENDS ${DICT})
endif()




set(LIBS ${ROOT_LIBRARIES} MathMore ) 

if (ANITA3_EVENTREADER) 
  add_definitions(-DANITA3_EVENTREADER)
  set(LIBS ${LIBS} AnitaEvent ) 
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
    link_directories(${AITA_UTIL_LIB_DIR})
    include_directories(${AITA_UTIL_LIB_DIR})
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


##We need to compile stuff into a library unless we want it to be compiled a
##million times over 
add_library(icemc_objs STATIC ${SOURCES} ${DICT}.cxx) 

macro(add_binary bin) 
  add_executable(${bin} ${bin}.cc)
  target_link_libraries(${bin} icemc_objs ${LIBS})
endmacro() 

add_binary(icemc) 
add_binary(testTrigger) 
add_binary(testSettings) 
add_binary(testEAS) 
add_binary(testThermalNoise)
add_binary(testInputAfterAntenna)
