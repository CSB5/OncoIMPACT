################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/lib/annotator_calculator.cpp \
../src/lib/annotator_input_output.cpp \
../src/lib/driver_genes.cpp \
../src/lib/explained_genes.cpp \
../src/lib/impact_scores.cpp \
../src/lib/intput.cpp \
../src/lib/merge_and_trim.cpp \
../src/lib/parameters.cpp \
../src/lib/phenotype_genes.cpp \
../src/lib/results.cpp \
../src/lib/sampling.cpp \
../src/lib/utilities.cpp 

OBJS += \
./src/lib/annotator_calculator.o \
./src/lib/annotator_input_output.o \
./src/lib/driver_genes.o \
./src/lib/explained_genes.o \
./src/lib/impact_scores.o \
./src/lib/intput.o \
./src/lib/merge_and_trim.o \
./src/lib/parameters.o \
./src/lib/phenotype_genes.o \
./src/lib/results.o \
./src/lib/sampling.o \
./src/lib/utilities.o 

CPP_DEPS += \
./src/lib/annotator_calculator.d \
./src/lib/annotator_input_output.d \
./src/lib/driver_genes.d \
./src/lib/explained_genes.d \
./src/lib/impact_scores.d \
./src/lib/intput.d \
./src/lib/merge_and_trim.d \
./src/lib/parameters.d \
./src/lib/phenotype_genes.d \
./src/lib/results.d \
./src/lib/sampling.d \
./src/lib/utilities.d 


# Each subdirectory must supply rules for building sources it contributes
src/lib/%.o: ../src/lib/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -std=c++0x -O3 -g3 -Wall -c -fmessage-length=0 -fopenmp -static-libgcc -static-libstdc++ -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


