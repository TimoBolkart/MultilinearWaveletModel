#
# Just finds the relevant include directory and libraries used to develop
# with glew

FIND_PATH(GLEW_INCLUDE_DIR glew.h
  /usr/local/include
  /usr/include
  /opt/local/include
)

FIND_LIBRARY(GLEW_LIBRARY
  NAMES
  glew
  PATHS
   /usr/local/lib
   /usr/lib
   /opt/local/lib
)

SET(GLEW_FOUND 0)
IF(GLEW_INCLUDE_DIR)
  IF(GLEW_LIBRARY)
    SET(GLEW_FOUND 1)
  ENDIF(GLEW_LIBRARY)
ENDIF(GLEW_INCLUDE_DIR)

IF(GLEW_FOUND)
  INCLUDE_DIRECTORIES(${GLEW_INCLUDE_DIR})
ELSE(GLEW_FOUND)
  MESSAGE("PROBLEM: Glew not found.")
ENDIF(GLEW_FOUND)

MARK_AS_ADVANCED(GLEW_INCLUDE_DIR 
  GLEW_LIBRARY)
