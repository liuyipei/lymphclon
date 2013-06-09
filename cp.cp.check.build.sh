rm -rf lymphclon
R CMD BATCH package.lymphclon.R
cp lymphclon_man/* lymphclon/man/
cp DESCRIPTION lymphclon/
R CMD check lymphclon
R CMD build lymphclon
