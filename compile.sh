# g++ -O3 -I ./  -lsuitesparseconfig -lumfpack -lamd -lcholmod -lcolamd -lcamd -lccolamd main.cpp -o solve 
g++ -std=gnu++17 -Wall -O3 -g -o solve main.cpp Simplex.cpp GS.cpp Data.cpp mpsReader.cpp Scaling.cpp -I /usr/include/suitesparse -lumfpack -lcholmod -lamd -lsuitesparseconfig
