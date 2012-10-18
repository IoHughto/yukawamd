#
#                           Makefile for MD tools
#                      ==============================
# 
#  This makefile requires a machine definition (.def) file. You can select one
#  of the following as default, or specify the MDEF macro on the make command
#  line.
#

#MDEF = AIX
 MDEF = BigRed
#MDEF = PGI
#MDEF = INTEL
#MDEF = MDGRAPE
#MDEF = Ranger_PGI


include $(MDEF).def

 FC = $(F95)
 FCFLAGS = $(F95OPT) $(TARCH) $(F95FLAGS)
 CCFLAGS = $(COPT) $(TARCH) $(CFLAGS)

 COMMON_OBJ = file_io.o file_stat.o
 NEW_FILE_RW = ftoc.o file_stat.o

 MAKEFILE = Makefile

dummy:
	@ echo -e "\nUsing definitions from $(MDEF).def\n"


#===================================================================================
#   Executables

xv8_xyz: xv8_xyz.o $(NEW_FILE_RW)
#$(CC) -o xv8_xyz xv8_xyz.o $(NEW_FILE_RW)
	$(FC) $(TARCH) $(FNOMAIN) -o xv8_xyz xv8_xyz.o $(NEW_FILE_RW)
	mv xv8_xyz /N/dc/scratch/jhughto/prog/bin/
#	make clean

vaf.3CP: vaf.3CP.o $(COMMON_OBJ)
#$(CC) -o vaf.3CP vaf.3CP.o $(COMMON_OBJ)
	$(FC) $(TARCH) $(FNOMAIN) -o vaf.3CP vaf.3CP.o $(COMMON_OBJ)
	mv vaf.3CP /N/dc/scratch/jhughto/prog/bin/
#	make clean

makebcc: makebcc.o $(COMMON_OBJ)
#$(CC) -o makebcc makebcc.o $(COMMON_OBJ)
	$(FC) $(TARCH) $(FNOMAIN) -o makebcc makebcc.o $(COMMON_OBJ)
	mv makebcc /N/dc/scratch/jhughto/prog/bin/
#	make clean

vaf.3CP.mpi: vaf.3CP.mpi.o $(COMMON_OBJ)
	mpif90 $(TARCH) $(FNOMAIN) -o vaf.3CP.mpi vaf.3CP.mpi.o $(COMMON_OBJ)
	mv vaf.3CP.mpi /N/dc/scratch/jhughto/prog/bin/
#	make clean

vaf.helium.mpi: vaf.helium.mpi.o $(COMMON_OBJ)
	mpif90 $(TARCH) $(FNOMAIN) -o vaf.helium.mpi vaf.helium.mpi.o $(COMMON_OBJ)
	mv vaf.helium.mpi /N/dc/scratch/jhughto/prog/bin/
#	make clean

vaf.3CP.mpi.2: vaf.3CP.mpi.2.o $(COMMON_OBJ)
	mpif90 $(TARCH) $(FNOMAIN) -o vaf.3CP.mpi.2 vaf.3CP.mpi.2.o $(COMMON_OBJ)
	mv vaf.3CP.mpi.2 /N/dc/scratch/jhughto/prog/bin/
#	make clean

vaf.3CP.mpi.double: vaf.3CP.mpi.double.o $(COMMON_OBJ)
	mpif90 $(TARCH) $(FNOMAIN) -o vaf.3CP.mpi.double vaf.3CP.mpi.double.o $(COMMON_OBJ)
	mv vaf.3CP.mpi.double /N/dc/scratch/jhughto/prog/bin/
#	make clean

solid.diffusion: solid.diffusion.o $(COMMON_OBJ)
#$(CC) -o solid.diffusion solid.diffusion.o $(COMMON_OBJ)
	$(FC) $(TARCH) $(FNOMAIN) -o solid.diffusion solid.diffusion.o $(COMMON_OBJ)
	mv solid.diffusion /N/dc/scratch/jhughto/prog/bin/
#	make clean

solid.diffusion.3CP: solid.diffusion.3CP.o $(COMMON_OBJ)
#$(CC) -o solid.diffusion.3CP solid.diffusion.3CP.o $(COMMON_OBJ)
	$(FC) $(TARCH) $(FNOMAIN) -o solid.diffusion.3CP solid.diffusion.3CP.o $(COMMON_OBJ)
	mv solid.diffusion.3CP /N/dc/scratch/jhughto/prog/bin/
#	make clean

solid.diffusion.3CP.full: solid.diffusion.3CP.full.o $(COMMON_OBJ)
#$(CC) -o solid.diffusion.3CP.full solid.diffusion.3CP.full.o $(COMMON_OBJ)
	$(FC) $(TARCH) $(FNOMAIN) -o solid.diffusion.3CP.full solid.diffusion.3CP.full.o $(COMMON_OBJ)
	mv solid.diffusion.3CP.full /N/dc/scratch/jhughto/prog/bin/
#	make clean

solid.diffusion.full: solid.diffusion.full.o $(COMMON_OBJ)
#$(CC) -o solid.diffusion.full solid.diffusion.full.o $(COMMON_OBJ)
	$(FC) $(TARCH) $(FNOMAIN) -o solid.diffusion.full solid.diffusion.full.o $(COMMON_OBJ)
	mv solid.diffusion.full /N/dc/scratch/jhughto/prog/bin/
#	make clean

solid.diffusion.stdev: solid.diffusion.stdev.o $(COMMON_OBJ)
#$(CC) -o solid.diffusion.stdev solid.diffusion.stdev.o $(COMMON_OBJ)
	$(FC) $(TARCH) $(FNOMAIN) -o solid.diffusion.stdev solid.diffusion.stdev.o $(COMMON_OBJ)
	mv solid.diffusion.stdev /N/dc/scratch/jhughto/prog/bin/
#	make clean

solid.diffusion.stdev.out: solid.diffusion.stdev.out.o $(COMMON_OBJ)
#$(CC) -o solid.diffusion.stdev.out solid.diffusion.stdev.out.o $(COMMON_OBJ)
	$(FC) $(TARCH) $(FNOMAIN) -o solid.diffusion.stdev.out solid.diffusion.stdev.out.o $(COMMON_OBJ)
	mv solid.diffusion.stdev.out /N/dc/scratch/jhughto/prog/bin/
#	make clean

diffusion.solid.2CP: diffusion.solid.2CP.o $(COMMON_OBJ)
#$(CC) -o diffusion.solid.2CP diffusion.solid.2CP.o $(COMMON_OBJ)
	$(FC) $(TARCH) $(FNOMAIN) -o diffusion.solid.2CP diffusion.solid.2CP.o $(COMMON_OBJ)
	mv diffusion.solid.2CP /N/dc/scratch/jhughto/prog/bin/
#	make clean

diffusion.solid: diffusion.solid.o $(COMMON_OBJ)
#$(CC) -o diffusion.solid diffusion.solid.o $(COMMON_OBJ)
	$(FC) $(TARCH) $(FNOMAIN) -o diffusion.solid diffusion.solid.o $(COMMON_OBJ)
	mv diffusion.solid /N/dc/scratch/jhughto/prog/bin/
#	make clean

s_q: s_q.o $(COMMON_OBJ)
#$(CC) -o s_q s_q.o $(COMMON_OBJ)
	$(FC) $(TARCH) $(FNOMAIN) -o s_q s_q.o $(COMMON_OBJ)
	mv s_q /N/dc/scratch/jhughto/prog/bin/
#	make clean

s_q.mpi: s_q.mpi.o $(COMMON_OBJ)
	mpif90 $(TARCH) $(FNOMAIN) -o s_q.mpi s_q.mpi.o $(COMMON_OBJ)
	mv s_q.mpi /N/dc/scratch/jhughto/prog/bin/
#	make clean

diffusion.vs.z: diffusion.vs.z.o $(COMMON_OBJ)
#$(CC) -o diffusion.vs.z diffusion.vs.z.o $(COMMON_OBJ)
	$(FC) $(TARCH) $(FNOMAIN) -o diffusion.vs.z diffusion.vs.z.o $(COMMON_OBJ)
	mv diffusion.vs.z /N/dc/scratch/jhughto/prog/bin/
#	make clean

nearest.neighbor: nearest.neighbor.o $(COMMON_OBJ)
#$(CC) -o nearest.neighbor nearest.neighbor.o $(COMMON_OBJ)
	$(FC) $(TARCH) $(FNOMAIN) -o nearest.neighbor nearest.neighbor.o $(COMMON_OBJ)
	mv nearest.neighbor /N/dc/scratch/jhughto/prog/bin/
#	make clean

nearest.neighbor.amorphous: nearest.neighbor.amorphous.o $(COMMON_OBJ)
#$(CC) -o nearest.neighbor.amorphous nearest.neighbor.amorphous.o $(COMMON_OBJ)
	$(FC) $(TARCH) $(FNOMAIN) -o nearest.neighbor.amorphous nearest.neighbor.amorphous.o $(COMMON_OBJ)
	mv nearest.neighbor.amorphous /N/dc/scratch/jhughto/prog/bin/
#	make clean

order: order.o $(COMMON_OBJ)
#$(CC) -o order order.o $(COMMON_OBJ)
	$(FC) $(TARCH) $(FNOMAIN) -o order order.o $(COMMON_OBJ)
	mv order /N/dc/scratch/jhughto/prog/bin/
#	make clean

order.stretch: order.stretch.o $(COMMON_OBJ)
#$(CC) -o order.stretch order.stretch.o $(COMMON_OBJ)
	$(FC) $(TARCH) $(FNOMAIN) -o order.stretch order.stretch.o $(COMMON_OBJ)
	mv order.stretch /N/dc/scratch/jhughto/prog/bin/
#	make clean

order.single: order.single.o $(COMMON_OBJ)
#$(CC) -o order.single order.single.o $(COMMON_OBJ)
	$(FC) $(TARCH) $(FNOMAIN) -o order.single order.single.o $(COMMON_OBJ)
	mv order.single /N/dc/scratch/jhughto/prog/bin/
#	make clean

nuc_statistics: nuc_statistics.o $(COMMON_OBJ)
#$(CC) -o nuc_statistics nuc_statistics.o $(COMMON_OBJ)
	$(FC) $(TARCH) $(FNOMAIN) -o nuc_statistics nuc_statistics.o $(COMMON_OBJ)
	mv nuc_statistics /N/dc/scratch/jhughto/prog/bin/
#	make clean

interface: interface.o $(COMMON_OBJ)
#$(CC) -o interface interface.o $(COMMON_OBJ)
	$(FC) $(TARCH) $(FNOMAIN) -o interface interface.o $(COMMON_OBJ)
	mv interface /N/dc/scratch/jhughto/prog/bin/
#	make clean

nucleation: nucleation.o $(COMMON_OBJ)
#$(CC) -o nucleation nucleation.o $(COMMON_OBJ)
	$(FC) $(TARCH) $(FNOMAIN) -o nucleation nucleation.o $(COMMON_OBJ)
	mv nucleation /N/dc/scratch/jhughto/prog/bin/
#	make clean

MCPdiff: MCPdiff.o $(COMMON_OBJ)
#$(CC) -o MCPdiff MCPdiff.o $(COMMON_OBJ)
	$(FC) $(TARCH) $(FNOMAIN) -o MCPdiff MCPdiff.o $(COMMON_OBJ)
	mv MCPdiff /N/dc/scratch/jhughto/prog/bin/
#	make clean

comp.bins: comp.bins.o $(COMMON_OBJ)
#$(CC) -o comp.bins comp.bins.o $(COMMON_OBJ)
	$(FC) $(TARCH) $(FNOMAIN) -o comp.bins comp.bins.o $(COMMON_OBJ)
	mv comp.bins /N/dc/scratch/jhughto/prog/bin/
#	make clean

phase.post: phase.post.o $(COMMON_OBJ)
#$(CC) -o phase.post phase.post.o $(COMMON_OBJ)
	$(FC) $(TARCH) $(FNOMAIN) -o phase.post phase.post.o $(COMMON_OBJ)
	mv phase.post /N/dc/scratch/jhughto/prog/bin/
#	make clean

fast.movers: fast.movers.o $(COMMON_OBJ)
#$(CC) -o fast.movers fast.movers.o $(COMMON_OBJ)
	$(FC) $(TARCH) $(FNOMAIN) -o fast.movers fast.movers.o $(COMMON_OBJ)
	mv fast.movers /N/dc/scratch/jhughto/prog/bin/
#	make clean

fast.movers.mpi: fast.movers.mpi.o $(COMMON_OBJ)
#$(CC) -o fast.movers.mpi fast.movers.mpi.o $(COMMON_OBJ)
	mpif90 $(TARCH) $(FNOMAIN) -o fast.movers.mpi fast.movers.mpi.o $(COMMON_OBJ)
	mv fast.movers.mpi /N/dc/scratch/jhughto/prog/bin/
#	make clean

shearmod.don: shearmod.don.o $(COMMON_OBJ)
#$(CC) -o shearmod.don shearmod.don.o $(COMMON_OBJ)
	$(FC) $(TARCH) $(FNOMAIN) -o shearmod.don shearmod.don.o $(COMMON_OBJ)
	mv shearmod.don /N/dc/scratch/jhughto/prog/bin/
#	make clean

shearmod.new: shearmod.new.o $(COMMON_OBJ)
#$(CC) -o shearmod.new shearmod.new.o $(COMMON_OBJ)
	$(FC) $(TARCH) $(FNOMAIN) -o shearmod.new shearmod.new.o $(COMMON_OBJ)
	mv shearmod.new /N/dc/scratch/jhughto/prog/bin/
#	make clean

shearmod.don.mpi: shearmod.don.mpi.o $(COMMON_OBJ)
	mpif90 $(TARCH) $(FNOMAIN) -o shearmod.don.mpi shearmod.don.mpi.o $(COMMON_OBJ)
	mv shearmod.don.mpi /N/dc/scratch/jhughto/prog/bin/
#	make clean

gofr: gofr.o $(COMMON_OBJ)
#$(CC) -o gofr gofr.o $(COMMON_OBJ)
	$(FC) $(TARCH) $(FNOMAIN) -o gofr gofr.o $(COMMON_OBJ)
	mv gofr /N/dc/scratch/jhughto/prog/bin/

gofr.mpi: gofr.mpi.o $(COMMON_OBJ)
	mpif90 $(TARCH) $(FNOMAIN) -o gofr.mpi gofr.mpi.o $(COMMON_OBJ)
	mv gofr.mpi /N/dc/scratch/jhughto/prog/bin/

gofr.mpi.full: gofr.mpi.full.o $(COMMON_OBJ)
	mpif90 $(TARCH) $(FNOMAIN) -o gofr.mpi.full gofr.mpi.full.o $(COMMON_OBJ)
	mv gofr.mpi.full /N/dc/scratch/jhughto/prog/bin/
#	make clean

gofr.mpi.stretch: gofr.mpi.stretch.o $(COMMON_OBJ)
	mpif90 $(TARCH) $(FNOMAIN) -o gofr.mpi.stretch gofr.mpi.stretch.o $(COMMON_OBJ)
	mv gofr.mpi.stretch /N/dc/scratch/jhughto/prog/bin/
#	make clean

3456to27648: 3456to27648.o $(COMMON_OBJ)
	$(FC) $(TARCH) $(FNOMAIN) -o 3456to27648 3456to27648.o $(COMMON_OBJ)
	mv 3456to27648 /N/dc/scratch/jhughto/prog/bin/

27648to55296: 27648to55296.o $(COMMON_OBJ)
	$(FC) $(TARCH) $(FNOMAIN) -o 27648to55296 27648to55296.o $(COMMON_OBJ)
	mv 27648to55296 /N/dc/scratch/jhughto/prog/bin/

target: $(MODOBJ) $(OBJ) $(EXTRAOBJ)
	$(FC) $(FLDFLAGS) -o $(TARGET) $(OBJ) $(MODOBJ) $(EXTRAOBJ) $(LIBS)

a2b: a2b.o $(COMMON_OBJ)
#$(CC) -o a2b a2b.o $(COMMON_OBJ)
	$(FC) $(TARCH) $(FNOMAIN) -o a2b a2b.o $(COMMON_OBJ)

composition.vs.z: composition.vs.z.o $(COMMON_OBJ)
#$(CC) -o composition.vs.z composition.vs.z.o $(COMMON_OBJ)
	$(FC) $(TARCH) $(FNOMAIN) -o composition.vs.z composition.vs.z.o $(COMMON_OBJ)
	mv composition.vs.z /N/dc/scratch/jhughto/prog/bin/
#	make clean

makeperfectbcc: makeperfectbcc.o $(COMMON_OBJ)
#$(CC) -o makeperfectbcc makeperfectbcc.o $(COMMON_OBJ)
	$(FC) $(TARCH) $(FNOMAIN) -o makeperfectbcc makeperfectbcc.o $(COMMON_OBJ)
	mv makeperfectbcc /N/dc/scratch/jhughto/prog/bin/
#	make clean

bench01: bench01.o
	$(FC) $(TARCH) -o bench01 bench01.o

ccat: ccat.o $(COMMON_OBJ)
	$(CC) $(TARCH) -o ccat ccat.o $(COMMON_OBJ) $(F95LIB)

clusteri: clusteri.o $(COMMON_OBJ)
	$(FC) $(TARCH) -o clusteri clusteri.o $(COMMON_OBJ)

cstat: cstat.o $(COMMON_OBJ)
#$(CC) -o cstat cstat.o $(COMMON_OBJ) $(F95LIB) -lm
	$(FC) $(TARCH) $(FNOMAIN) -o cstat cstat.o $(COMMON_OBJ)

diffconfig: diffconfig.o force_err.o $(COMMON_OBJ)
	$(FC) $(TARCH) -o diffconfig diffconfig.o force_err.o $(COMMON_OBJ)

dkb01: dkb01.o
	$(FC) $(TARCH) -o dkb01 dkb01.o

dkb02: dkb02.o $(COMMON_OBJ) randnums.o
	$(FC) $(TARCH) -o dkb02 dkb02.o $(COMMON_OBJ) randnums.o

formfactor: formfactor.o
	$(FC) $(TARCH) -o formfactor formfactor.o

hcalc: hcalc.o
	$(FC) $(TARCH) -o hcalc hcalc.o

jobstats: jobstats.o
	$(CC) $(TARCH) -o jobstats jobstats.o

mpi_test: mpi_test.o
	$(FC) $(TARCH) -o mpi_test mpi_test.o

pastelog: pastelog.o
	$(CC) $(TARCH) -o pastelog pastelog.o

pickout: pickout.o
	$(CC) $(TARCH) -o pickout pickout.o

rpl: rpl.o $(COMMON_OBJ)
#$(CC) -o rpl rpl.o $(COMMON_OBJ) $(F95LIB)
	$(FC) $(TARCH) $(FNOMAIN) -o rpl rpl.o $(COMMON_OBJ)

shft: shft.o $(COMMON_OBJ)
	$(CC) $(TARCH) $(FNOMAIN) -o shft shft.o $(COMMON_OBJ)

reduce: reduce.o $(COMMON_OBJ)
	$(FC) $(TARCH) $(FNOMAIN) -o reduce reduce.o $(COMMON_OBJ)

strata: strata.o $(COMMON_OBJ)
	$(FC) $(TARCH) -o strata strata.o $(COMMON_OBJ)


#===================================================================================
#   Main program files

a2b.o: a2b.c
	$(CC) -c $(CCFLAGS) a2b.c

bench01.o: bench01.f
	$(FC) -c $(FCFLAGS) bench01.f

cluster.o: cluster.f
	$(FC) -c $(FCFLAGS) cluster.f

clusteri.o: clusteri.f
	$(FC) -c $(FCFLAGS) clusteri.f

ccat.o: ccat.c
	$(CC) -c $(CCFLAGS) ccat.c

cstat.o: cstat.c
	$(CC) -c $(CCFLAGS) cstat.c

diffconfig.o: diffconfig.F
	$(FC) -c $(FCFLAGS) diffconfig.F $(GETARG)

dkb01.o: dkb01.c
	$(CC) -c $(CCFLAGS) dkb01.c

dkb02.o: dkb02.F
#$(FC) -c $(FCFLAGS) -WF,-DGETARG dkb02.F
	$(FC) -c $(FCFLAGS) -DGETARG dkb02.F

formfactor.o: formfactor.f
	$(FC) -c $(FCFLAGS) formfactor.f

hcalc.o: hcalc.f
	$(FC) -c $(FCFLAGS) hcalc.f

jobstats.o: jobstats.c
	$(CC) -c $(CCFLAGS) jobstats.c

mpi_test.o: mpi_test.f
	$(FC) -c $(FCFLAGS) mpi_test.f

pastelog.o: pastelog.c
	$(CC) -c $(CCFLAGS) pastelog.c

pickout.o: pickout.c
	$(CC) -c $(CCFLAGS) pickout.c

rpl.o: rpl.c
	$(CC) -c $(CCFLAGS) rpl.c

shft.o: shft.c
	$(CC) -c $(CCFLAGS) shft.c

strata.o: strata.F
	$(FC) -c $(FCFLAGS) -WF,-DGETARG strata.F
#$(FC) -c $(FCFLAGS) -DGETARG strata.F


#-----------------------------------------------------------------
#   Subprogram files

force_err.o: force_err.f
	$(FC) -c $(FCFLAGS) force_err.f

ftoc.o: ftoc.f
	$(FC) -c $(FCFLAGS) ftoc.f

file_io.o: file_io.f
	$(FC) -c $(FCFLAGS) file_io.f

file_stat.o: file_stat.c
	$(CC) -c $(CCFLAGS) file_stat.c

vaf.3CP.mpi.o: vaf.3CP.mpi.c
	$(MPCC) -c $(CCFLAGS) vaf.3CP.mpi.c

vaf.helium.mpi.o: vaf.helium.mpi.c
	$(MPCC) -c $(CCFLAGS) vaf.helium.mpi.c

vaf.3CP.mpi.2.o: vaf.3CP.mpi.2.c
	$(MPCC) -c $(CCFLAGS) vaf.3CP.mpi.2.c

vaf.3CP.mpi.double.o: vaf.3CP.mpi.double.c
	$(MPCC) -c $(CCFLAGS) vaf.3CP.mpi.double.c

s_q.mpi.o: s_q.mpi.c
	$(MPCC) -c $(CCFLAGS) s_q.mpi.c

shearmod.don.mpi.o: shearmod.don.mpi.c
	$(MPCC) -c $(CCFLAGS) shearmod.don.mpi.c

gofr.mpi.o: gofr.mpi.c
	$(MPCC) -c $(CCFLAGS) gofr.mpi.c

gofr.mpi.full.o: gofr.mpi.full.c
	$(MPCC) -c $(CCFLAGS) gofr.mpi.full.c

gofr.mpi.stretch.o: gofr.mpi.stretch.c
	$(MPCC) -c $(CCFLAGS) gofr.mpi.stretch.c

fast.movers.mpi.o: fast.movers.mpi.c
	$(MPCC) -c $(CCFLAGS) fast.movers.mpi.c

randnums.o: randnums.f
	$(FC) -c $(FCFLAGS) randnums.f



#===================================================================================
clean:
	rm -f *.o

realclean:
	rm -f *.o a2b ccat cluster cstat diffconfig rpl shft
