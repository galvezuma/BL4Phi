################################################################################
# SGR
################################################################################

-include ../makefile.init

RM := rm -rf
CC := icl
CFLAGS := /Qmic -O3 -c 
LFLAGS := /Qmic -o
LIBS := -lpthread -lrt

# All of the sources participating in the build are defined here
-include sources.mk
-include src/subdir.mk
-include subdir.mk
-include objects.mk

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
endif

-include ../makefile.defs

# Add inputs and outputs from these tool invocations to the build variables 

# All Target


# Tool invocations
all: 
	@echo 'Building'
	pwd
	$(CC) $(CFLAGS) Farrar.c $(LIBS) -o ../Release/Farrar.o
	$(CC) $(CFLAGS) GeneralFunctions.c $(LIBS) -o ../Release/GeneralFunctions.o
	$(CC) $(CFLAGS) ManageDatabase.c $(LIBS) -o ../Release/ManageDatabase.o
	$(CC) $(CFLAGS) MultipleQuery.c $(LIBS) -o ../Release/MultipleQuery.o
	$(CC) $(CFLAGS) SingleQuery.c $(LIBS) -o ../Release/SingleQuery.o
	$(CC) $(CFLAGS) BLVector.c $(LIBS) -o ../Release/BLVector.o
	@echo 'Linking'
	$(CC) $(LFLAGS) ../Release/BLVector $(LIBS) ../Release/BLVector.o ../Release/SingleQuery.o ../Release/MultipleQuery.o ../Release/ManageDatabase.o ../Release/GeneralFunctions.o ../Release/Farrar.o

# Other Targets
clean:
	-$(RM) $(OBJS)$(C_DEPS)$(EXECUTABLES) BLVector
	-@echo 'Cleaned'
	pwd

.PHONY: all clean dependents
.SECONDARY:

-include ../makefile.targets
