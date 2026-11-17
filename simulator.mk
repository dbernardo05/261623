COMPILER = g++
CCFLAGS = -std=c++11 -g -fopenmp
COMMON_HEADERS = ../neuron_output.h ../output_options.h ../run_options.h ../run_output.h ../snapshot.h ../membrane_noise.h

macgregor: main.o layer.o neuron.o simulation.o
	$(COMPILER) $(CCFLAGS) -o $@ $^

main.o: main.cpp layer.h neuron.h ../run_options.h ../output_options.h
	$(COMPILER) $(CCFLAGS) -c $<
layer.o: layer.cpp layer.h neuron.h ../neuron_output.h ../output_options.h ../run_options.h
	$(COMPILER) $(CCFLAGS) -c $<
neuron.o: neuron.cpp neuron.h ../membrane_noise.h
	$(COMPILER) $(CCFLAGS) -c $<
simulation.o: ../simulation.cpp layer.h neuron.h $(COMMON_HEADERS)
	$(COMPILER) $(CCFLAGS) -I. -c $< -o $@

.PHONY: reference clean
reference: macgregor-reference
macgregor-reference: ../reference/build.sh ../reference/recording.patch ../simulation.cpp $(COMMON_HEADERS) ../simulator.mk
	sh ../reference/build.sh $(notdir $(CURDIR)) "$(COMPILER)" $(CCFLAGS)

clean:
	$(RM) *.o macgregor macgregor-reference
