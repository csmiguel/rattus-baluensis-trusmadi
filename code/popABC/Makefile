DIR = src
CC = gcc
CFLAGS = -std=gnu89 -Wall -O2
LIBS = -lm

OBJS_POP = \
 $(patsubst %.c,%.o,$(DIR)/interface.c) $(patsubst %.c,%.o,$(DIR)/mylib.c) $(patsubst %.c,%.o,$(DIR)/pop_abc.c) \
 $(patsubst %.c,%.o,$(DIR)/pop_convertabc2.c) $(patsubst %.c,%.o,$(DIR)/pop_convertabc3.c) $(patsubst %.c,%.o,$(DIR)/pop_convertabc4.c) \
 $(patsubst %.c,%.o,$(DIR)/pop_convertabc.c) $(patsubst %.c,%.o,$(DIR)/pop_firstpass.c) $(patsubst %.c,%.o,$(DIR)/pop_genetictree.c) \
 $(patsubst %.c,%.o,$(DIR)/pop_joindata.c) $(patsubst %.c,%.o,$(DIR)/pop_makepop.c) $(patsubst %.c,%.o,$(DIR)/pop_makeprior.c) \
 $(patsubst %.c,%.o,$(DIR)/pop_makestats.c) $(patsubst %.c,%.o,$(DIR)/pop_maketarget.c) $(patsubst %.c,%.o,$(DIR)/pop_samplepriors.c) \
 $(patsubst %.c,%.o,$(DIR)/pop_summstats.c)

OBJS_REJ = \
 $(patsubst %.c,%.o,$(DIR)/firstpass.c) $(patsubst %.c,%.o,$(DIR)/mylib.c)

OBJS_SHU = \
 $(patsubst %.c,%.o,$(DIR)/mylib.c) $(patsubst %.c,%.o,$(DIR)/turnit.c)

OBJS_SIM = \
 $(patsubst %.c,%.o,$(DIR)/abc.c) $(patsubst %.c,%.o,$(DIR)/genetictree.c) $(patsubst %.c,%.o,$(DIR)/mylib.c) \
 $(patsubst %.c,%.o,$(DIR)/samplepriors.c) $(patsubst %.c,%.o,$(DIR)/summstats.c)

OBJS_SUM = \
 $(patsubst %.c,%.o,$(DIR)/maketarget.c) $(patsubst %.c,%.o,$(DIR)/mylib.c) $(patsubst %.c,%.o,$(DIR)/summstats.c)

all: popabc rejection shuffle simulate summdata

popabc: $(OBJS_POP)
	$(CC) $(CFLAGS) $(OBJS_POP) -o $@ $(LIBS)

rejection: $(OBJS_REJ)
	$(CC) $(CFLAGS) $(OBJS_REJ) -o $@ $(LIBS)

shuffle: $(OBJS_SHU)
	$(CC) $(CFLAGS) $(OBJS_SHU) -o $@ $(LIBS)

simulate: $(OBJS_SIM)
	$(CC) $(CFLAGS) $(OBJS_SIM) -o $@ $(LIBS)

summdata: $(OBJS_SUM)
	$(CC) $(CFLAGS) $(OBJS_SUM) -o $@ $(LIBS)

clean:
	rm -f $(wildcard $(DIR)/*.o) popabc rejection shuffle simulate summdata
