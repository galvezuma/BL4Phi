SRCS = Farrar.c GeneralFunctions.c ManageDatabase.c MultipleQuery.c SingleQuery.c BLPhi.c
OBJS = $(SRCS:.c=.o)

#
# Please, note the usage of /Qmic for Xeon-Phi. Remove when compiling for
# other architectures. In addition, the option /Q is used when using the 
# Intel compiler in Windows.
#

RM := rm -rf
CC := icl
CFLAGS := /Qmic -O3 -c 
LFLAGS := /Qmic -o
LIBS := -lpthread -lrt
MAIN = BLVector

all: $(MAIN)

$(MAIN): $(OBJS) 
	@echo 'Linking'
	$(CC) $(LFLAGS) $(MAIN) $(LIBS) $(OBJS)
	

.c.o:
	$(CC) $(CFLAGS) $^ $(LIBS) -o $@

.PHONY: clean

clean:
	$(RM) $(OBJS) $(MAIN)
	
