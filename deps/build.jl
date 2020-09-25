println("building lib/f16_fortran")
cd("../lib/f16_fortran")
run(`mkdir -p bin`)
run(`make`)
