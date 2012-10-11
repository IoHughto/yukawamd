 MAKEFILE = Makefile

#================================================================================
#

# GNU compilers
#F77        = g77
#CC         = gcc

# Select a machine definition file to define compilers and options.
#MDEF       = AIX
#MDEF       = INTEL
#MDEF       = PGI
 MDEF       = BigRed
#MDEF       = Quarry

include ../mdef/$(MDEF).def

 FCFLAGS    = $(F95OPT) $(TARCH) $(F95FLAGS) $(ATP)
 FCVERSION  = $(F95VERSION)
 CCFLAGS    = $(TARCH)

 FLDFLAGS   = $(F95OPT) $(TARCH)

 LDFLAGS    =
#LDFLAGS    = -L/usr/local/lib

 MODS       = md_types.mod md_globals.mod md_comm.mod
 MODOBJ     = md_types.o md_globals.o
 OBJ        = md.o input1.o initconfig.o randnums.o strain.o int_control.o md_utils.o \
    gen_utils.o rw_xva.o rw_xvfa.o rw_xvb.o rw_xvfb.o file_stat.o perf.o \
    dwallclock00.o $(ATPOBJ)
 INTOBJGPC  = int_nn.o int_ion_pure.o int_ion_mix.o
 INTEGRATOR = vel_verlet_2.o

dummy:
	@ echo -e "\nUsing definitions from $(MDEF).def\n"

#===============================================================================
#      Executables

md: dummy
	$(F95) $(F95VERSION)
	@ echo -e ""
	cp mpif_stubs.h mpif.h
	make -f $(MAKEFILE) target "TARGET=md" "FC=$(F95)" "CC=$(CC)" \
              "MODS=$(MODULES)" \
              "MODOBJ=$(MODOBJ) md_comm_ser.o" "EXTRAOBJ=$(INTOBJGPC) $(INTEGRATOR)" \
              "LIBS="

md_omp: dummy
	cp mpif_stubs.h mpif.h
	make -f $(MAKEFILE) target "TARGET=md_omp" "FC=$(OMPF95)" "CC=$(OMPCC)" \
              "MODOBJ=$(MODOBJ) md_comm_ser.o" "EXTRAOBJ=$(INTOBJGPC) $(INTEGRATOR)" \
              "LIBS="

md_mpi: dummy
	rm -f mpif.h
	make -f $(MAKEFILE) target "TARGET=md_mpi" "FC=$(MPF95)" "CC=$(MPCC)" \
              "MODOBJ=$(MODOBJ) md_comm_mpi.o" "EXTRAOBJ=$(INTOBJGPC) $(INTEGRATOR)" \
              "LIBS="

md_mpi_omp: dummy
	rm -f mpif.h
	make -f $(MAKEFILE) target "TARGET=md_mpi_omp" "FC=$(MPOMPF95)" "CC=$(MPOMPCC)" \
              "MODOBJ=$(MODOBJ) md_comm_mpi.o" "EXTRAOBJ=$(INTOBJGPC) $(INTEGRATOR)" \
              "LIBS="


#===============================================================================
#  Utilities and other programs.

mcat: mcat.o
	$(F95) $(FLDFLAGS) -o mcat mcat.o


#===============================================================================
#     Each make of a program sets several macros, and then recursively calls
#     make. The target of this recursive call is "target":

target: $(MODOBJ) $(OBJ) $(EXTRAOBJ)
	$(FC) $(FLDFLAGS) -o $(TARGET) $(OBJ) $(MODOBJ) $(EXTRAOBJ) $(LIBS)


#===============================================================================
#     Object files

md_types.o: md_types.f
	$(FC) -c $(FCFLAGS) md_types.f

md_globals.o: md_types.mod md_globals.f
	$(FC) -c $(FCFLAGS) md_globals.f

md_comm_ser.o: md_globals.mod md_comm_ser.f
	$(FC) -c $(FCFLAGS) md_comm_ser.f

md_comm_mpi.o: md_globals.mod md_comm_mpi.f
	$(FC) -c $(FCFLAGS) md_comm_mpi.f

md.o: $(MODS) md.F perf.h
	$(FC) -c $(FCFLAGS) md.F

input1.o: $(MODS) input1.F
	$(FC) -c $(FCFLAGS) input1.F

initconfig.o: $(MODS) initconfig.f
	$(FC) -c $(FCFLAGS) initconfig.f

vel_verlet_1.o: $(MODS) vel_verlet_1.f perf.h
	$(FC) -c $(FCFLAGS) vel_verlet_1.f

vel_verlet_2.o: $(MODS) vel_verlet_2.f perf.h
	$(FC) -c $(FCFLAGS) vel_verlet_2.f

vel_verlet_3.o: $(MODS) vel_verlet_3.f perf.h
	$(FC) -c $(FCFLAGS) vel_verlet_3.f

strain.o: $(MODS) strain.f perf.h
	$(FC) -c $(FCFLAGS) strain.f

int_control.o: $(MODS) int_control.f perf.h
	$(FC) -c $(FCFLAGS) int_control.f

int_nn.o: $(MODS) int_nn.f perf.h
	$(FC) -c $(FCFLAGS) int_nn.f

int_ion_pure.o: $(MODS) int_ion_pure.f perf.h
	$(FC) -c $(FCFLAGS) int_ion_pure.f

int_ion_mix.o: $(MODS) int_ion_mix.f perf.h
	$(FC) -c $(FCFLAGS) int_ion_mix.f

md_utils.o: md_utils.f
	$(FC) -c $(FCFLAGS) md_utils.f

gen_utils.o: $(MODS) gen_utils.f
	$(FC) -c $(FCFLAGS) gen_utils.f

rw_xva.o: $(MODS) rw_xva.f
	$(FC) -c $(FCFLAGS) rw_xva.f

rw_xvfa.o: $(MODS) rw_xvfa.f
	$(FC) -c $(FCFLAGS) rw_xvfa.f

rw_xvb.o: $(MODS) rw_xvb.f
	$(FC) -c $(FCFLAGS) rw_xvb.f

rw_xvfb.o: $(MODS) rw_xvfb.f
	$(FC) -c $(FCFLAGS) rw_xvfb.f

randnums.o: randnums.f
	$(FC) -c $(FCFLAGS) randnums.f

g.o: $(MODS) g.f
	$(FC) -c $(FCFLAGS) g.f

s.o: $(MODS) s.f
	$(FC) -c $(FCFLAGS) s.f

perf.o: perf.f perf.h
	$(FC) -c $(FCFLAGS) perf.f

file_stat.o: file_stat.c
	$(CC) -c $(CCFLAGS) file_stat.c

dwallclock00.o: dwallclock00.c
	$(CC) -c $(CCFLAGS) dwallclock00.c

ATP_diags.o: ATP_diags.f
	$(FC) -c $(FCFLAGS) ATP_diags.f

mcat.o: mcat.f
	$(F95) -c $(FCFLAGS) mcat.f

clean:
	rm -f *.o *.mod

realclean:
	rm -f *.o *.mod \
              md md_omp md_mpi md_mpi_omp \
