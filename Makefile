PLATFORM=icc

include Makefile.in.$(PLATFORM)

.PHONY:	all
all: $(OBJS)
	$(LD) $(OBJS) -o interpolate $(LDFLAGS) $(LIBS) $(OMP_CFLAGS)

%.o: %.c
	$(CC) -c $(CFLAGS) $(OPTFLAGS) $(CPPFLAGS) $(OMP_CFLAGS)$<

.PHONY:	clean
clean:
	rm -f *.o
