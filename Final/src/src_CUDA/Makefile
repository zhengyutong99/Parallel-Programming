CC = g++
NVCC = nvcc
CCFLAG = -Wall -std=c++17 -g -Wno-unused-result 
NVCCFLAG= -std=c++17 -gencode=arch=compute_86,code=sm_86 -Xcompiler '-fPIC'
OUTPUT = bestdock_gemdock_CUDA
# OUTPUT = bestdock_gemdock

OBJ = gemdock2.o kernel.o 

all: $(OUTPUT)

$(OUTPUT): $(OBJ)
	$(NVCC) $(OBJ) -o $(OUTPUT)

%.o: %.cpp
	$(CC) $(CCFLAG) -c $< -o $@

%.o: %.cu
	$(NVCC) $(NVCCFLAG) -c $< -o $@

clean:
	rm -f *.o $(OUTPUT)