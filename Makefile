.PHONY: docs serve build deploy clean install

install:
	pip install mkdocs mkdocs-material

PORT ?= 8000

serve:
	mkdocs serve -a 127.0.0.1:$(PORT)

build:
	mkdocs build

deploy:
	mkdocs gh-deploy

clean:
	rm -rf site/
