install:
	python setup.py install

test:
	nosetests -v

clean:
	rm -rf build/
	rm -rf dist/
	rm -rf gwass.egg-info/	 	
