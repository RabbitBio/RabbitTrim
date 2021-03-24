bin/ktrim: src/ktrim.cpp src/common.h src/util.h src/param_handler.h src/pe_handler.h src/se_handler.h src/io/*.h src/find_seed.h src/io/*.cpp
	@echo Build Ktrim
	cd src/io; g++ -std=c++11 -c FastxStream.cpp -o FastxStream.o -g -O3;cd ..; cd ..
	cd src/io; g++ -std=c++11 -c Formater.cpp -o Formater.o -g -O3; cd ..; cd ..
	cd src; g++ -std=c++11 -c ktrim.cpp -o ktrim.o -g -mavx512f -mavx512bw -O3; cd ..
	cd src; g++ -std=c++11 -fopenmp -lz -mavx512f -mavx512bw -I ./io/ -O3 ktrim.o ./io/FastxStream.o ./io/Formater.o -o ../bin/ktrim; cd ..

install: bin/ktrim	# requires root
	@echo Install Ktrim for all users
	@cp bin/ktrim /usr/local/bin

clean:
	rm -f bin/ktrim

