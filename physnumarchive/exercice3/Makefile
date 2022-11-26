CC = g++
CFLAGS = -std=c++0x -Wall -g
EXEC_NAME_STD = Exercice3
INCLUDES =
LIBS =
OBJ_FILES = Exercice3_mysolution.o
EXEC_SUFFIX = .exe

# add suffix to exec
EXEC_NAME_SUFFIX = $(addsuffix $(EXEC_SUFFIX),$(EXEC_NAME))
EXEC_NAME_STD_SUFFIX = $(addsuffix $(EXEC_SUFFIX),$(EXEC_NAME_STD))

clean :
	rm $(EXEC_NAME_SUFFIX) $(EXEC_NAME_STD_SUFFIX)  $(OBJ_FILES) *.out
 
solution : $(OBJ_FILES)
	$(CC) -o $(EXEC_NAME_SUFFIX) $(OBJ_FILES) $(LIBS)

%.o: %.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

