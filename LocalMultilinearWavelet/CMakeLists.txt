SET( LocalMultilinearWaveletFiles
	BSGWFactory.cpp
	FaceData.cpp
	GlobalImageMatch.cpp
	FileLoader.cpp
	MMProjectionCostFunction_Local.cpp
	MultilinearModel.cpp
	MultilinearModelHandler.cpp
	NearestNeighborAssistant.cpp
	Profile.cpp
	StereoFace.cpp
	WaveletShapeMultiLinearOptimizer.cpp
	BSGWFactory.h
	BSplineGridWavelet.h
	FaceData.h
	GlobalImageMatch.h
	FileLoader.h
	MMProjectionCostFunction_Local.h
	MultilinearModel.h
	MultilinearModelHandler.h
	NearestNeighborAssistant.h
	Profile.h
	StereoFace.h
	WaveletShapeFitter.h
	WaveletShapeMultiLinearOptimizer.h
)

INCLUDE_DIRECTORIES(${ITK_INCLUDES}) 

INCLUDE_DIRECTORIES(${ANN_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${CLAPACK_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${OpenCV_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${GLUT_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${GLEW_INCLUDE_DIR})

ADD_LIBRARY(LocalMultilinearWaveletLib ${LocalMultilinearWaveletFiles})
TARGET_LINK_LIBRARIES(LocalMultilinearWaveletLib ${ITK_LIBRARIES} ${ANN_LIBRARIES} ${CLAPACK_LIBRARIES} ${OpenCV_LIBS} ${GLUT_LIBRARY} ${GLEW_LIBRARY})

ADD_EXECUTABLE(LocalMultilinearWavelet LocalMultilinearWavelet.cpp)
TARGET_LINK_LIBRARIES(LocalMultilinearWavelet LocalMultilinearWaveletLib)



	