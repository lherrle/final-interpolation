PLATFORM=icc

include Makefile.in.$(PLATFORM)

.PHONY:	all
all: $(OBJS)
	$(CC) $(OBJS) -o interpolate $(LDFLAGS) $(LIBS)

%.o: %.c
	$(CC) -c $(CFLAGS) $(OPTFLAGS) $(CPPFLAGS) $(OMP_CFLAGS)$<

.PHONY:	clean
clean:
	rm -f *.o
