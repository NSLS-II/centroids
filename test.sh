#!/bin/bash -

set -o nounset -e                                 # Treat unset variables as an error

while [ 1 ]; do
	py.test -s -k "test_find_photons"
done
