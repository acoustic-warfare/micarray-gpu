CC = gcc
CFLAGS = -v -Wall -Werror -Os
CUDA_LIBS = -L/usr/local/cuda/lib64

RM = rm

SRC = src
OUT = build

.PHONY: all

all: cubridge.so

bridge.o:
	$(CC) $(CFLAGS) -I /usr/local/cuda/include -c -o $(OUT)/bridge.o $(SRC)/bridge.c -lstdc++

beam.o:
	nvcc --compiler-options '-fPIC' -c -o $(OUT)/beam.o $(SRC)/beam.cu

common.o:
	$(CC) $(CFLAGS) -c -o $(OUT)/common.o $(SRC)/common.c

cubridge.so: bridge.o beam.o common.o
	$(CC) $(CFLAGS) -shared -o $(OUT)/cubridge.so $(OUT)/bridge.o $(OUT)/beam.o $(OUT)/common.o $(CUDA_LIBS) -lcudart -lstdc++ -fPIC


.PHONY: clean

clean:
	$(RM) $(OUT)/*.o $(OUT)/*.so