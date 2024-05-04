CXX      ?= g++
CXXFLAGS ?= -std=c++20
CPPFLAGS ?= -O3 -Wall -Wno-conversion-null -Wno-deprecated-declarations 



# Executable name
EXEC = main

# Source files
SRCS = main.cpp
# Object files
OBJS = $(SRCS:.cpp=.o)

all: $(EXEC)

%.o: %.cpp 
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@

$(EXEC): $(OBJS)
	$(CXX) $(OBJS) $(CXXFLAGS) $(LDFLAGS) -o $@

clean:
	$(RM) *.o $(EXEC)

distclean: clean