test-build:
	docker build . -f binder/Dockerfile --build-arg HOME=/root