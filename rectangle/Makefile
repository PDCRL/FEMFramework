#Open terminal in the same folder and run make. (type make and enter)
# change target names as you like. it is the executable. you can run executable in terminal by typing ./Oceane 

TARGET = Oceane

LIBS = -lm -lnlopt

#don't change anything below.
SRC = $(wildcard ../src/*.cpp)
SRC += main.cpp
OBJ = $(patsubst src/%.cpp, build/%.o, $(SRC))

OFILES = $(patsubst src/%.cpp,build/%.o,$(SRC))
CXXFLAGS =  -I../include -I/usr/local/include
LFLAGS = -L/usr/local/lib -lnlopt -lm 

$(TARGET):$(SRC)
	g++ $^   $(CXXFLAGS) $(LFLAGS) -o $@
 
					
clean :
	rm -r build $(TARGET)

	



