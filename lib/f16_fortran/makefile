FC = gfortran
FFLAGS = -Wall -fdefault-real-8 -fdefault-double-8 -Wline-truncation -shared -fPIC
SWITCH = -O3

SRC_PATH = ./src
BIN_PATH = ./bin

SRC = $(wildcard $(SRC_PATH)/*.f)

all: so

so:
	$(FC) $(SWITCH) $(FFLAGS) $(SRC) -o $(BIN_PATH)/f16.so

clean:
	rm -f bin/*
	rm -f build/*
