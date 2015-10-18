################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/annotator.cpp \
../src/database.cpp \
../src/discovery.cpp \
../src/oncoIMPACT.cpp 

OBJS += \
./src/annotator.o \
./src/database.o \
./src/discovery.o \
./src/oncoIMPACT.o 

CPP_DEPS += \
./src/annotator.d \
./src/database.d \
./src/discovery.d \
./src/oncoIMPACT.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -std=c++0x -O3 -g3 -Wall -c -fmessage-length=0 -fopenmp -static-libgcc -static-libstdc++ -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


