# This software is Copyright (C) 2004-2008  Bristol University
# and is released under the GNU General Public License version 2.

all:	makecpp makerun

# For now don't optimise python or
# make Mathematica dir
# makepy makema

makepy:	./py/
	cd ./py/ ; make

makecpp:	./cpp/
	cd ./cpp/ ; jam

makerun:	./config-run
	cd ./config-run/ ; jam

makema:	./ma/
	cd ./ma/ ; make

clean:	
	cd ./cpp/ ; jam clean
	cd ./config-run/ ; jam clean
	cd ./py/ ; make clean
	cd ./ma/ ; make clean

