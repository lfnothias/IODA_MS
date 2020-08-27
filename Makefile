test-build:
	docker build . -f binder/Dockerfile --build-arg HOME=/root

build-standalone:
	docker build . -t ioda-standalone -f docker/Dockerfile

run_standalone:
	docker run -p 9000:8888 --rm -it ioda-standalone /bin/bash

run_standalone_notebook:
	docker run -p 9000:8888 -it ioda-standalone

test-unit:
	act -P ubuntu-latest=nektos/act-environments-ubuntu:18.04 -b

