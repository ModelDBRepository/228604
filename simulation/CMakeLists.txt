cmake_minimum_required (VERSION 2.6)

set(PROJECTNAME sn)
set(SPONT_EXENAME sn_spont)
set(SIG_EXENAME sn_sig)
set(FPT_EXENAME sn_fpt)
set(QDA_EXENAME qifda)

project (${PROJECTNAME})

add_custom_command (
	OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/autoparams.h ${CMAKE_CURRENT_BINARY_DIR}/autoparams.cpp 
	COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/param_gen.py ${CMAKE_CURRENT_SOURCE_DIR}/parameters
	DEPENDS parameters
)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)

execute_process(WORKING_DIRECTORY "@CMAKE_CURRENT_SOURCE_DIR@" 
                COMMAND git status --porcelain
                OUTPUT_VARIABLE GIT_STAT_OUT)
if ( NOT ${GIT_STAT_OUT} STREQUAL "")
    set(GIT_SHA1 "${GIT_SHA1}_dirty")
endif()

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/git_info.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/git_info.cpp" @ONLY)

find_library(LIB_BOOST_PO boost_program_options)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_definitions(-DBOOST_TEST_DYN_LINK) 

# include git rev in name of executable
set(SPONT_EXENAME ${SPONT_EXENAME}_${GIT_SHA1})
set(SIG_EXENAME ${SIG_EXENAME}_${GIT_SHA1})
set(FPT_EXENAME ${FPT_EXENAME}_${GIT_SHA1})
set(QDA_EXENAME ${QDA_EXENAME}_${GIT_SHA1})

add_executable(${SPONT_EXENAME} sn_spont.cpp neurophys/fourier_transform.cpp neurophys/spike_train.cpp neurophys/signal.cpp neurophys/array_funcs.cpp neurophys/histogram_funcs.cpp ${CMAKE_CURRENT_BINARY_DIR}/autoparams.cpp)
target_link_libraries(${SPONT_EXENAME} boost_program_options gsl fftw3 gslcblas)
add_executable(${SIG_EXENAME} sn_sig.cpp neurophys/fourier_transform.cpp neurophys/spike_train.cpp neurophys/signal.cpp neurophys/array_funcs.cpp neurophys/histogram_funcs.cpp ${CMAKE_CURRENT_BINARY_DIR}/autoparams.cpp)
target_link_libraries(${SIG_EXENAME} boost_program_options gsl fftw3 gslcblas)
add_executable(${FPT_EXENAME} fpt_gsl_sn.cpp ${CMAKE_CURRENT_BINARY_DIR}/autoparams.cpp)
target_link_libraries(${FPT_EXENAME} boost_program_options gsl fftw3 gslcblas)
add_executable(${QDA_EXENAME} qif_diffap.cpp ${CMAKE_CURRENT_BINARY_DIR}/autoparams.cpp)
target_link_libraries(${QDA_EXENAME} boost_program_options gsl fftw3 gslcblas)

add_custom_command(
    TARGET ${SPONT_EXENAME}_${GIT_SHA1} POST_BUILD
    COMMAND ln -sf ${SPONT_EXENAME}_${GIT_SHA1} ${SPONT_EXENAME}
)
add_custom_command(
    TARGET ${SIG_EXENAME}_${GIT_SHA1} POST_BUILD
    COMMAND ln -sf ${SIG_EXENAME}_${GIT_SHA1} ${SIG_EXENAME}
)
