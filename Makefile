all: visionTestPNG

makeVisionTest : visionTestPNG.cpp
	clang++ $< -o $@
