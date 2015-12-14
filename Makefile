PLATFORM=icc

include Makefile.in.$(PLATFORM)

.PHONY:	all
all: $(OBJS)
	$(LD) $(OBJS) -o interpolate $(LDFLAGS) $(LIBS)

%.o: %.c
	$(CC) -c $(CFLAGS) $(OPTFLAGS) $(CPPFLAGS) $<

.PHONY:	clean
clean:
	rm -f *.o
