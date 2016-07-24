# Makefile for SUPOLAR_PS Tests
# Copied from ...su/main

include $(CWPROOT)/src/Makefile.config


D = $L/libcwp.a $L/libpar.a $L/libsu.a


LFLAGS= $(PRELFLAGS) -L$L -L/usr/include -lgsl -lgslcblas -lsu -lpar -lcwp -lm -fopenmp $(POSTLFLAGS)

PROGS =			\
supolar_PS
#        sustalta
#	$B/sualford	\
#	$B/sueipofi	\
#	$B/suhrot	\
#	$B/sultt	\
#	$B/supofilt	\
#	$B/supolar


INSTALL	:	$(PROGS)
	@-rm -f INSTALL
	@touch $@


$(PROGS):	$(CTARGET) $D 
	-$(CC) $(CFLAGS) $(@F).c $(LFLAGS) -o $@
	@$(MCHMODLINE)
	@echo $(@F) installed in $B

remake	:
	-rm -f $(PROGS) INSTALL
	$(MAKE) 
	
clean::
	rm -f a.out junk* JUNK* core
