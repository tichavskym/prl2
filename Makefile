PROC=1
STEPS=2


.PHONY: build run clean

all: build

build:
	mpic++ --prefix /usr/local/share/OpenMPI -o life life.cpp

run: build
	mpirun --prefix /usr/local/share/OpenMPI --use-hwthread-cpus -np $(PROC) life grid.txt $(STEPS)

clean:
	rm life
#generate:
#	dd if=/dev/random bs=1 count=$(NUMBERS) of=numbers 2> /dev/null

#run: build
#	# TODO use-hwthread-cpus is my addition to the command from assignment
#	mpirun --prefix /usr/local/share/OpenMPI --use-hwthread-cpus -np $(PROC) pms
#
#clean:
#	rm numbers pms