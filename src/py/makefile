# This software is Copyright (C) 2004-2008  Bristol University
# and is released under the GNU General Public License version 2.

all:	compile check lint test

compile:
	python2.3 ./ImportAllOptimized.py -O .

check:
	pychecker *.py

lint:
	pylint *.py

test:

clean:	
	-rm degreeHij
	-rm degreeWi
	-rm degreeKi
	-rm *~
	-rm *.pyo
	-rm *.pyc
	-rm ../tmp/*_example*.db
	-rm ../tmp/*_test*.db
	-rm /tmp/*9.db
	-rm /tmp/*8.db
	-rm /tmp/*7.db
	-rm /tmp/*6.db
	-rm /tmp/*5.db
	-rm /tmp/*4.db
	-rm /tmp/*3.db
	-rm /tmp/*2.db
	-rm /tmp/*1.db
	-rm /tmp/*0.db

