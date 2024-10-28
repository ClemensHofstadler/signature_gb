.PHONY: install install-ahocorasick build clean

install: build install-ahocorasick
	sage -pip install --upgrade .

install-ahocorasick:
	AHOCORASICK_BYTES=yes \
	sage -pip install ./pyahocorasick-master

build:
	sage -python setup.py build_ext --inplace
	
clean:
	rm -rf __pycache__