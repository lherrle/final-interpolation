PLATFORM=icc

include Makefile.in.$(PLATFORM)

.PHONY:	all
all: $(OBJS)
	$(CC) $(OBJS) -o interpolate

%.o: %.c
	$(CC) -c $(CFLAGS) $(OPTFLAGS) $(CPPFLAGS) $<

.PHONY:	clean
clean:
	rm -f *.o
