PLATFORM=icc

include Makefile.in.$(PLATFORM)

.PHONY:	all
all:	$(OBJS)

%.o: %.c
$(CC) -c $(CFLAGS) $(OPTFLAGS) $(CPPFLAGS) $<