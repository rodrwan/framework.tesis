mexcmd = 'mex -outdir bin';
mexcmd = [mexcmd ' -O'];
mexcmd = [mexcmd ' CXXOPTIMFLAGS="-O3 -DNDEBUG"'];
mexcmd = [mexcmd ' LDOPTIMFLAGS="-O3"'];
mexcmd = [mexcmd ' CXXFLAGS="\$CXXFLAGS -Wall"'];
mexcmd = [mexcmd ' LDFLAGS="\$LDFLAGS -Wall"'];

eval([mexcmd ' features/features.cc -o features']);