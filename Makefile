.PHONY: help build test lint clean docs shell batch

help:
	@echo "Available commands:"
	@echo "  make test          - Run the CLI test"
	@echo "  make build         - Build the Docker image"
	@echo "  make build-test     - Run the CLI test in Docker"
	@echo "  make shell         - Open a shell inside the Docker container"
	@echo "  make docs          - Generate HTML documentation (if Sphinx is configured)"
	@echo "  make clean         - Remove build artifacts and caches"

test:
	uv run --dev pre-commit install && uv run --dev pre-commit run --all
	uv run --dev coverage run -m pytest --junitxml=report.xml tests
	uv run --dev coverage xml -o coverage/cobertura-coverage.xml
	uv run --dev coverage report -m

build:
	docker build -t afdb-toolkit .

build-test:
	docker run --rm afdb-toolkit test


shell:
	docker run --rm -it afdb-toolkit /bin/bash

docs:
	uv run --group docs sphinx-build docs build/html

clean:
	rm -rf build dist *.egg-info .pytest_cache .mypy_cache .ruff_cache coverage .coverage report.xml coverage.xml coverage-html.xml
