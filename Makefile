
INC_DIR=/opt/local/include

# boost_program_options-xgcc42-mt-x64-1_70

BOOST_VERSION=1_70
BOOST_INC_ROOT=/opt/local/include/boost-$(BOOST_VERSION)
BOOST_LIB_DIR=/opt/local/lib
BOOST_SUFFIX=xgcc42-mt-x64-$(BOOST_VERSION)
BOOST_PROGRAM_OPTIONS_LIB=boost_program_options-$(BOOST_SUFFIX)

CXXFLAGS += -I$(BOOST_INC_ROOT) -I$(INC_DIR)
LFLAGS   += -L$(BOOST_LIB_DIR) -l$(BOOST_PROGRAM_OPTIONS_LIB) -lpng

all: visionTestPNG

visionTestPNG : visionTestPNG.cpp
	clang++ $(CXXFLAGS) $(LFLAGS) $< -o $@

testImage : visionTestPNG
	./visionTestPNG /tmp/a_test.png
