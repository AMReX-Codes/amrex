# Makefile created by mkmf 19.3.0


MK_TEMPLATE = template.mk
include $(MK_TEMPLATE)


.DEFAULT:
	-echo $@ does not exist.
all: main3d.MPI.gnu.ex
Adv_3d.o: Src_3d/Adv_3d.f90 compute_flux_3d.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	Src_3d/Adv_3d.f90
Prob.o: ../Exec/SingleVortex/Prob.f90
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../Exec/SingleVortex/Prob.f90
amr_data_mod.o: ./amr_data_mod.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./amr_data_mod.F90
averagedown_mod.o: ./averagedown_mod.F90 amr_data_mod.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./averagedown_mod.F90
bc_mod.o: ./bc_mod.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./bc_mod.F90
compute_dt_mod.o: ./compute_dt_mod.F90 my_amr_mod.o face_velocity_3d.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./compute_dt_mod.F90
compute_flux_3d.o: Src_3d/compute_flux_3d.f90 slope_3d.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	Src_3d/compute_flux_3d.f90
evolve_mod.o: ./evolve_mod.F90 my_amr_mod.o amr_data_mod.o compute_dt_mod.o plotfile_mod.o averagedown_mod.o face_velocity_3d.o Adv_3d.o fillpatch_mod.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./evolve_mod.F90
face_velocity_3d.o: ../Exec/SingleVortex/face_velocity_3d.f90
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../Exec/SingleVortex/face_velocity_3d.f90
fillpatch_mod.o: ./fillpatch_mod.F90 amr_data_mod.o bc_mod.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./fillpatch_mod.F90
fmain.o: ./fmain.F90 my_amr_mod.o initdata.o evolve_mod.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./fmain.F90
initdata.o: ./initdata.F90 my_amr_mod.o plotfile_mod.o averagedown_mod.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./initdata.F90
my_amr_mod.o: ./my_amr_mod.F90 amr_data_mod.o bc_mod.o Prob.o fillpatch_mod.o tagging_mod.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./my_amr_mod.F90
plotfile_mod.o: ./plotfile_mod.F90 my_amr_mod.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./plotfile_mod.F90
slope_3d.o: Src_3d/slope_3d.f90
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	Src_3d/slope_3d.f90
tagging_mod.o: ./tagging_mod.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./tagging_mod.F90
./compute_flux_3d.f90: Src_3d/compute_flux_3d.f90
	cp Src_3d/compute_flux_3d.f90 .
./slope_3d.f90: Src_3d/slope_3d.f90
	cp Src_3d/slope_3d.f90 .
./Adv_3d.f90: Src_3d/Adv_3d.f90
	cp Src_3d/Adv_3d.f90 .
./face_velocity_3d.f90: ../Exec/SingleVortex/face_velocity_3d.f90
	cp ../Exec/SingleVortex/face_velocity_3d.f90 .
./Prob.f90: ../Exec/SingleVortex/Prob.f90
	cp ../Exec/SingleVortex/Prob.f90 .
SRC = Src_3d/compute_flux_3d.f90 ../Exec/SingleVortex/face_velocity_3d.f90 ./fillpatch_mod.F90 ./tagging_mod.F90 ./plotfile_mod.F90 ./fmain.F90 ./initdata.F90 ./compute_dt_mod.F90 Src_3d/slope_3d.f90 ./amr_data_mod.F90 ./evolve_mod.F90 Src_3d/Adv_3d.f90 ../Exec/SingleVortex/Prob.f90 ./averagedown_mod.F90 ./bc_mod.F90 ./my_amr_mod.F90
OBJ = compute_flux_3d.o face_velocity_3d.o fillpatch_mod.o tagging_mod.o plotfile_mod.o fmain.o initdata.o compute_dt_mod.o slope_3d.o amr_data_mod.o evolve_mod.o Adv_3d.o Prob.o averagedown_mod.o bc_mod.o my_amr_mod.o
OFF = Src_3d/compute_flux_3d.f90 Src_3d/slope_3d.f90 Src_3d/Adv_3d.f90 ../Exec/SingleVortex/face_velocity_3d.f90 ../Exec/SingleVortex/Prob.f90
clean: neat
	-rm -f .main3d.MPI.gnu.ex.cppdefs $(OBJ) main3d.MPI.gnu.ex *.mod

neat:
	-rm -f $(TMPFILES)

localize: $(OFF)
	cp $(OFF) .

TAGS: $(SRC)
	etags $(SRC)

tags: $(SRC)
	ctags $(SRC)

main3d.MPI.gnu.ex: $(OBJ)
	$(LD) $(OBJ) -o main3d.MPI.gnu.ex  $(LDFLAGS)

