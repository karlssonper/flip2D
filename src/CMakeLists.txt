PROJECT(FLIP2D_SRC)

SET(SOURCE flip2D grid particles sdf pressure pcg multigrid gaussSeidel jacobi)

ADD_LIBRARY(flip2D SHARED ${SOURCE})

TARGET_LINK_LIBRARIES(flip2D)

INSTALL(TARGETS flip2D DESTINATION lib)