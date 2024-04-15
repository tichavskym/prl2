PROC=1
STEPS=4


.PHONY: build run clean

all: build

build:
	mpic++ --prefix /usr/local/share/OpenMPI -o life life.cpp

run: build
	mpirun --prefix /usr/local/share/OpenMPI --use-hwthread-cpus -np $(PROC) life grid.txt $(STEPS)

clean:
	rm life
