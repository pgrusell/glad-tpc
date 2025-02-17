##############################################################################
#   Copyright (C) 2018 GSI Helmholtzzentrum für Schwerionenforschung GmbH    #
#   Copyright (C) 2018-2025 Members of the R3B Collaboration                 #
#                                                                            #
#             This software is distributed under the terms of the            #
#                 GNU Lesser General Public Licence (GPL) version 3,         #
#                    copied verbatim in the file "LICENSE".                  #
#                                                                            #
# In applying this license GSI does not waive the privileges and immunities  #
# granted to it by virtue of its status as an Intergovernmental Organization #
# or submit itself to any jurisdiction.                                      #
##############################################################################

# Create a library called "libR3BGTPCDATA" which includes the source files given 
# in the array. The extension is already found. Any number of sources could be 
# listed here.

set(SRCS
    R3BGTPCPoint.cxx
    R3BGTPCProjPoint.cxx
    R3BGTPCMappedData.cxx
    R3BGTPCCalData.cxx
    R3BGTPCHitData.cxx
    R3BGTPCHitClusterData.cxx
    R3BGTPCTrackData.cxx)

# fill list of header files from list of source files
# by exchanging the file extension
change_file_extension(*.cxx *.h HEADERS "${SRCS}")

add_library_with_dictionary(
    LIBNAME
    R3BGTPCData
    LINKDEF
    R3BGTPCDataLinkDef.h
    HEADERS
    ${HEADERS}
    SRCS
    ${SRCS}
    INCLUDEDIRS
    ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDENCIES
    R3BTracking)

target_compile_options(R3BGTPCData PRIVATE -Wno-overloaded-virtual)
