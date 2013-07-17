rm -rf lymphclon
R CMD BATCH package.lymphclon.R
cp lymphclon_man/* lymphclon/man/
mkdir lymphclon/man/
cp DESCRIPTION lymphclon/
cp NAMESPACE lymphclon/
R CMD check lymphclon
R CMD build lymphclon
rm *.Rout
