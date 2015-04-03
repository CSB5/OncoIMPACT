################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/driver_genes.cpp \
../src/explained_genes.cpp \
../src/intput.cpp \
../src/oncoIMPACT.cpp \
../src/phenotype_genes.cpp \
../src/sampling.cpp \
../src/tester.cpp \
../src/utilities.cpp 

OBJS += \
./src/driver_genes.o \
./src/explained_genes.o \
./src/intput.o \
./src/oncoIMPACT.o \
./src/phenotype_genes.o \
./src/sampling.o \
./src/tester.o \
./src/utilities.o 

CPP_DEPS += \
./src/driver_genes.d \
./src/explained_genes.d \
./src/intput.d \
./src/oncoIMPACT.d \
./src/phenotype_genes.d \
./src/sampling.d \
./src/tester.d \
./src/utilities.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -std=c++11 -D__NO_INLINE__ -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


