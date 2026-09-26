PROJ_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
EXTENSION_NAME=duckvep
USE_UNSTABLE_C_API=0
TARGET_DUCKDB_VERSION=v1.2.0
DUCKDB_HEADER_VERSION=v1.5.3

ifeq ($(DUCKDB_PLATFORM),windows_amd64_mingw)
override GEN=
override VCPKG_TOOLCHAIN_PATH=
override VCPKG_TARGET_TRIPLET=
override VCPKG_HOST_TRIPLET=
endif
ifeq ($(DUCKDB_PLATFORM),windows_amd64_rtools)
override GEN=
override VCPKG_TOOLCHAIN_PATH=
override VCPKG_TARGET_TRIPLET=
override VCPKG_HOST_TRIPLET=
endif

include extension-ci-tools/makefiles/c_api_extensions/base.Makefile
include extension-ci-tools/makefiles/c_api_extensions/c_cpp.Makefile

.PHONY: all test test_debug test_release test-extension-symbols readme
all: configure release
configure: venv platform extension_version
platform: venv
extension_version: venv platform
	@$(VERSION_COMMAND)
configure/extension_version.txt: venv platform
check_configure: configure
build_extension_with_metadata_debug build_extension_with_metadata_release: extension_version
debug: build_extension_library_debug build_extension_with_metadata_debug
release: build_extension_library_release build_extension_with_metadata_release
# CI builds in its container and then runs test_release on the host against those
# artifacts, so the test targets must not rebuild. Local `make test` builds first.
test: debug
	$(MAKE) test_debug
test_debug: test_extension_debug
test_release: test_extension_release test-extension-symbols
test-extension-symbols:
	@set -e; file=build/release/duckvep.duckdb_extension; test -f "$$file"; \
		if test "$$(uname -s)" = Darwin; then \
			actual=$$(nm -gU "$$file" | awk '$$2 ~ /^[TDB]$$/ {sub(/^_/, "", $$3); print $$3}' | sort -u); \
		else \
			actual=$$(nm -D --defined-only "$$file" | awk '$$2 ~ /^[TDB]$$/ && $$3 !~ /^(_init|_fini)$$/ {print $$3}' | sort -u); \
			nm -D -u "$$file" | awk '$$NF ~ /^duckdb_/ {print "Unexpected DuckDB API import: " $$0; bad=1} END {exit bad}'; \
		fi; \
		test "$$actual" = duckvep_init_c_api || { printf 'Unexpected exports: %s\n' "$$actual"; exit 1; }
readme: release
	Rscript scripts/render-readme.R

.PHONY: site site-check
site:
	Rscript scripts/build-site.R

site-check: site
	Rscript scripts/check-site-links.R
