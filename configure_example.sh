# /usr/bin.bash

./configure -prefix=$SU2_RUN/.. --with-MPI=mpicxx --with-Metis-lib=$SU2_HOME/metis-4.0.3/ --with-Metis-include=$SU2_HOME/metis-4.0.3/Lib/

# ./configure --with-MPI=mpicxx --with-Metis-lib=/usr/local/metis-4.0.3 --with-Metis-include=/usr/local/metis-4.0.3/Lib --prefix=${SU2_HOME}

#./configure -prefix=$SU2_RUN/.. --with-Metis-lib=$SU2_HOME/metis-4.0.3/ --with-Metis-include=$SU2_HOME/metis-4.0.3/Lib/ 

# ./configure -prefix=$SU2_RUN/.. --with-MPI=mpicxx --with-Metis-lib=$SU2_HOME/metis-4.0.3/ --with-Metis-include=$SU2_HOME/metis-4.0.3/Lib/ --disable-CFD
