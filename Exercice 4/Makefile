CC = g++
CFLAGS = -std=c++0x -Wall -g
EXEC_NAME_APO = Exercice4_2022_Apollo_solution
EXEC_NAME_ATM = Exercice4_2022_atmos_solution 
INCLUDES =
LIBS =
OBJ_FILES_APO = Exercice4_2022_Apollo_solution.o
OBJ_FILES_ATM = Exercice4_2022_atmos_solution.o
EXEC_SUFFIX = .exe

# add suffix to exec
EXEC_NAME_APO_SUFFIX = $(addsuffix $(EXEC_SUFFIX),$(EXEC_NAME_APO))
EXEC_NAME_ATM_SUFFIX = $(addsuffix $(EXEC_SUFFIX),$(EXEC_NAME_ATM))

all : apollo, atmos

clean :
	rm $(EXEC_NAME_SUFFIX)  $(OBJ_FILES) *.out
 
apollo : $(OBJ_FILES_APO)
	$(CC) -o $(EXEC_NAME_APO_SUFFIX) $(OBJ_FILES_APO) $(LIBS)

atmos : $(OBJ_FILES_ATM)
	$(CC) -o $(EXEC_NAME_ATM_SUFFIX) $(OBJ_FILES_ATM) $(LIBS)

%.o: %.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

