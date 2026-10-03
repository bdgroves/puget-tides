# PUGET-TIDES - make (fetch + build + run), make build, make run,
# make normals (30 years of hourly water levels, a few minutes), make clean
FC     = gfortran
FFLAGS = -O2 -Wall -Wno-character-truncation

.PHONY: all fetch build run normals clean

all: fetch build run

fetch:
	python3 fetch_tides.py || [ $$? -lt 8 ]

build: puget-tides

puget-tides: TIDE_PHYS.f90 PUGET-TIDES.f90
	$(FC) $(FFLAGS) -o $@ $^

normals-bin: TIDE_PHYS.f90 NORMALS.f90
	$(FC) $(FFLAGS) -o normals $^

normals: normals-bin
	python3 fetch_tides.py --history
	./normals || [ $$? -lt 8 ]
	rm -rf history

run: puget-tides
	./puget-tides || [ $$? -lt 8 ]
	@cat puget-tides-report.txt

clean:
	rm -f puget-tides normals *.mod *.tmp water_level.csv pressure.csv monthly_mean.csv noaa_hilo.csv fetch_status.csv
