test-build:
	docker build . -f binder/Dockerfile --build-arg HOME=/root

test-unit:
	act -P ubuntu-latest=nektos/act-environments-ubuntu:18.04 -b

