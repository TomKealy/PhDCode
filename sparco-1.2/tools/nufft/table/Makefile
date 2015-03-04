# Makefile
# one can reduce time by compiling the initialization routine

suf_i686 = mexglx
suf_ppc = mexmac
suf_x86_64 = mexa64
suf_i386 = mexmaci
suf_ = suf_i686; # default is linux
suf = $(suf_$(Arch))
dir = ../../mex/v7/
def = def,table.h
mexloc = /usr/local/matlab/bin/mex
mexloc = mex

goal:	\
	$(dir)interp1_table_mex.$(suf) \
	$(dir)interp2_table_mex.$(suf) \
	$(dir)interp3_table_mex.$(suf) \
	$(dir)interp1_table_adj_mex.$(suf) \
	$(dir)interp2_table_adj_mex.$(suf) \
	$(dir)interp3_table_adj_mex.$(suf)
	@echo working on $(suf)

aflag_x86_64 = -fPIC
aflag = $(aflag_$(Arch))

# this -D_GNU_SOURCE flag came from looking at "mex -v" output...
CFLAGS = -D_GNU_SOURCE \
	-ansi -Wall -Wshadow -Wpointer-arith -W -Werror $(aflag)
mex = $(mexloc) CFLAGS='$(CFLAGS)'

i1c = interp1_table_mex.c
$(dir)interp1_table_mex.$(suf):	$(i1c) $(def)
	$(mex) -output $@ $(i1c)
	chmod 755 $@

i2c = interp2_table_mex.c
$(dir)interp2_table_mex.$(suf):	$(i2c) $(def)
	$(mex) -output $@ $(i2c)
	chmod 755 $@

i3c = interp3_table_mex.c
$(dir)interp3_table_mex.$(suf):	$(i3c) $(def)
	$(mex) -output $@ $(i3c)
	chmod 755 $@

a1c = interp1_table_adj_mex.c
$(dir)interp1_table_adj_mex.$(suf):	$(a1c) $(def)
	$(mex) -output $@ $(a1c)
	chmod 755 $@

a2c = interp2_table_adj_mex.c
$(dir)interp2_table_adj_mex.$(suf):	$(a2c) $(def)
	$(mex) -output $@ $(a2c)
	chmod 755 $@

a3c = interp3_table_adj_mex.c
$(dir)interp3_table_adj_mex.$(suf):	$(a3c) $(def)
	$(mex) -output $@ $(a3c)
	chmod 755 $@


reset:
	/bin/rm -i $(dir)interp*table*mex.$(suf)
