cmake_minimum_required(VERSION 3.12)

project(Scratchapixel CXX C)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

if(WIN32)
	set(PYTHON_EXECUTABLE "python")
else()
	set(PYTHON_EXECUTABLE "python3")
endif()

if((NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/vendor/.bootstrap.json") OR ("${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt" IS_NEWER_THAN "${CMAKE_CURRENT_SOURCE_DIR}/vendor/.bootstrap.json"))
	execute_process(
		COMMAND ${PYTHON_EXECUTABLE} "bootstrap.py"
		WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
		RESULT_VARIABLE bootstrap_result
	)
	if (NOT (${bootstrap_result} EQUAL 0))
		message(FATAL_ERROR "Cannot run bootstrap.py")
	endif()
endif()

if(WIN32)
	add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif()

include_directories(.)

# add samples
add_subdirectory(MinimalRayTracer/RayBox)
add_subdirectory(MinimalRayTracer/SimpleShapes)
add_subdirectory(Geometry/Points_Vectors_Normals)




