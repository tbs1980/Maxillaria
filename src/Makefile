MAXILLARIA_BUILD_TARGET=gcc.linux-g10
#MAXILLARIA_BUILD_TARGET=suncc.linux-g10
include Makefile.$(MAXILLARIA_BUILD_TARGET)

default:maxillaria.exe libmaxillaria.a

LIB_MAXILLARIA_OBJS=nested_sampler.o live_point.o mt19937.o nested_diagnostic.o priors.o utils.o post_equal_weights.o post_stats.o explore.o

MAXILLARIA_OBJS=maxillaria.o $(LIB_MAXILLARIA_OBJS)

libmaxillaria.a:$(LIB_MAXILLARIA_OBJS)
	$(AR) $@ $^
maxillaria.exe:$(MAXILLARIA_OBJS)
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)
%.o:%.c
	$(CC) $(CFLAGS) -c $<

clean:
	rm *.o *.exe *.a
	
distclean:
	rm *.o *.exe *.txt *.dat

