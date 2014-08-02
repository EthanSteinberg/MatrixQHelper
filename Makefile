ifeq (1, 0)
	CXX = x86_64-w64-mingw32-g++
	JNI_INCLUDE = /usr/lib/jvm/java-8-oracle/include/win32
	LIB_FOLDER = lib/win32
	LIB_PREFIX = 
	LIB_EXTENSION = .dll
	EXE_EXTENSION = .exe
	LIB_STATIC = -static
else
	CXX = g++
	JNI_INCLUDE = /usr/lib/jvm/java-8-oracle/include/linux
	LIB_FOLDER = lib/linux
	LIB_EXTENSION = .so
	LIB_PREFIX = lib
	EXE_EXTENSION = 
	LIB_STATIC = 
endif


SRCS = src/main.cpp src/jni_impl.cpp
DEPS = $(SRCS:.cpp=.d)
OBJECTS = $(SRCS:.cpp=.o)

CXXFLAGS = -std=c++11 -Wall -fPIC -march=native -ffast-math -O2
CPPFLAGS := -I /usr/lib/jvm/java-8-oracle/include -I ${JNI_INCLUDE} -I ./include

LFLAGS = -L ${LIB_FOLDER}

LIBS := -Wl,-Bstatic ${LIB_STATIC} -static-libgcc -static-libstdc++ 

EXECUTABLE := program${EXE_EXTENSION}
LIBRARY := ${LIB_PREFIX}JNIMatrixQHelper${LIB_EXTENSION}

default: ${LIBRARY} ${EXECUTABLE} matrixqhelperlib/JNIMatrixQHelper.class

-include $(DEPS)

%.d : %.cpp
	$(CXX) ${CPPFLAGS} $(CXXFLAGS) -MF"$@" -MM -MT"$@" -MT"$(<:.cpp=.o)" "$<"

${LIBRARY}: $(OBJECTS)
	$(CXX) -shared ${OBJECTS} ${LFLAGS} ${LIBS} -o $@ 

${EXECUTABLE}: $(OBJECTS)
	$(CXX) ${OBJECTS} ${LFLAGS} ${LIBS} -o $@

%.class: %.java
	javac  -source 1.7 -target 1.7 $<

clean:
	-rm -f $(OBJECTS)
	-rm -f $(DEPS)
	-rm -f $(LIBRARY)
	-rm -f $(EXECUTABLE)
	
