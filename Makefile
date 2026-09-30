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

# Scale-runner aggregator guards: synthetic receipts, no DuckVEP run.
.PHONY: test_scale_aggregate
test_scale_aggregate:
	scripts/duckvep_scale_selftest_aggregate.sh

.PHONY: all test test_debug test_release test_fault_injection build_fault_injection test_haplotype_contract test-extension-symbols test_mane_grch37 test_mane_grch37_receipt readme build_asan test_release_asan test_properties test_properties_sanitized
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
test_release: test_extension_release test-extension-symbols test-sql-lambda-syntax

# Fault injection: an AddressSanitizer + LeakSanitizer build whose allocation
# budget can fail the Nth allocation (duckvep_fault_arm). Every allocation made
# by model load and annotation is failed in turn; each failure must be a clean
# error that publishes nothing and leaks nothing. Linux/GCC-compatible hosts.
build_fault_injection: check_configure
	cmake $(CMAKE_VERSION_PARAMS) -DCMAKE_BUILD_TYPE=Debug \
		"-DCMAKE_C_FLAGS=-fsanitize=address -fno-omit-frame-pointer -O1 -g -DDUCKVEP_FAULT_INJECTION=1" \
		-DCMAKE_SHARED_LINKER_FLAGS=-fsanitize=address -S $(PROJ_DIR) -B cmake_build/fault
	cmake --build cmake_build/fault --parallel
	mkdir -p build/fault
	$(PYTHON_VENV_BIN) extension-ci-tools/scripts/append_extension_metadata.py \
		-l cmake_build/fault/$(EXTENSION_LIB_FILENAME) -o build/fault/$(EXTENSION_FILENAME) \
		-n $(EXTENSION_NAME) -dv $(TARGET_DUCKDB_VERSION) \
		-evf configure/extension_version.txt -pf configure/platform.txt
test_fault_injection: build_fault_injection
	python3 test/scripts/test_fault_injection.py

# Pure policy fixtures run offline. The full-release receipt needs the external
# Parquet output and is invoked explicitly after building a staged release.
test_mane_grch37:
	Rscript test/scripts/test_mane_grch37_policy.R
	@if test -n "$(MANE_GRCH37_OUTPUT)"; then $(MAKE) test_mane_grch37_receipt MANE_GRCH37_OUTPUT="$(MANE_GRCH37_OUTPUT)"; else echo 'Full-release receipt needs MANE_GRCH37_OUTPUT'; fi

test_haplotype_contract:
	Rscript --vanilla test/scripts/check_haplotype_csq_map.R test/data/haplotype/csq_so_map_v1.tsv
	Rscript --vanilla test/scripts/check_haplotype_goldens.R
	Rscript --vanilla test/scripts/generate_haplotype_models.R vertical same_codon frame startstop nmd
	Rscript --vanilla test/scripts/check_same_codon_goldens.R
	Rscript --vanilla test/scripts/check_frame_goldens.R
	Rscript --vanilla test/scripts/check_startstop_goldens.R
	Rscript --vanilla test/scripts/check_nmd_goldens.R
	Rscript --vanilla test/scripts/test_haplotype_accounting.R

# Native property suite (the host-neutral kernel, greatest + theft) under AddressSanitizer and
# UndefinedBehaviorSanitizer. Offline. Any sanitizer report aborts (no recovery, leak checking on) and
# any failed property is a nonzero exit, so the target fails on either. PROPERTY_CC and PROPERTY_FLAGS
# may be overridden (for example PROPERTY_CC=clang); DUCKVEP_PROPERTY_TRIALS/SEED tune the suite.
PROPERTY_CC ?= cc
PROPERTY_BUILD := $(PROJ_DIR)build/properties
PROPERTY_SANITIZE := -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer
PROPERTY_FLAGS ?= -std=gnu11 -O1 -g
PROPERTY_SOURCES := $(wildcard $(PROJ_DIR)test/duckvep/property/*.c) \
	$(wildcard $(PROJ_DIR)test/duckvep/vendor/theft/src/*.c) $(wildcard $(PROJ_DIR)src/kernel/src/*.c)
# Host-neutral C properties over the kernel views (theft + greatest, no DuckDB),
# always under ASan and UBSan. DUCKVEP_PROP_TRIALS and DUCKVEP_PROP_SEED make a
# run larger or reproducible; read_sources.pl checks every test is registered.
test_properties:
	perl $(PROJ_DIR)test/duckvep/property/read_sources.pl >/dev/null
	mkdir -p $(PROPERTY_BUILD)
	$(PROPERTY_CC) $(PROPERTY_FLAGS) $(PROPERTY_SANITIZE) \
		-I$(PROJ_DIR)src/kernel/include -I$(PROJ_DIR)src/kernel/src -I$(PROJ_DIR)test/duckvep/property \
		-I$(PROJ_DIR)test/duckvep/vendor/greatest -I$(PROJ_DIR)test/duckvep/vendor/theft/inc \
		-I$(PROJ_DIR)test/duckvep/vendor/theft/src \
		$(PROPERTY_SOURCES) -lm -o $(PROPERTY_BUILD)/duckvep_properties
	ASAN_OPTIONS=detect_leaks=1:abort_on_error=0:halt_on_error=1:strict_string_checks=1:detect_stack_use_after_return=1 \
	UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1 \
		$(PROPERTY_BUILD)/duckvep_properties

test_properties_sanitized: test_properties

test_mane_grch37_receipt:
	@test -n "$(MANE_GRCH37_OUTPUT)" || { echo 'Set MANE_GRCH37_OUTPUT to the full-release output directory' >&2; exit 1; }
	Rscript test/scripts/test_mane_grch37_receipt.R "$(MANE_GRCH37_OUTPUT)" benchmarks/data/mane_grch37_receipts.csv

# DuckDB 2.0 rejects single-arrow SQL lambdas by default.
.PHONY: test-sql-lambda-syntax
test-sql-lambda-syntax:
	python3 test/scripts/check_sql_lambdas.py
test-extension-symbols:
	@set -e; file=build/release/duckvep.duckdb_extension; test -f "$$file"; \
		case "$$(uname -s)" in MINGW*|MSYS*|CYGWIN*|Windows_NT) windows=1 ;; *) windows=0 ;; esac; \
		if test "$$(uname -s)" = Darwin; then \
			actual=$$(nm -gU "$$file" | awk '$$2 ~ /^[TDB]$$/ {sub(/^_/, "", $$3); print $$3}' | sort -u); \
		elif test "$$windows" = 1; then \
			actual=$$(objdump -p "$$file" | awk '/\[Ordinal\/Name Pointer\] Table/ {t=1; next} t && /^\t\[ *[0-9]+\]/ {print $$NF; next} t && NF == 0 {t=0}' | sort -u); \
		else \
			actual=$$(nm -D --defined-only "$$file" | awk '$$2 ~ /^[TDB]$$/ && $$3 !~ /^(_init|_fini)$$/ {print $$3}' | sort -u); \
			nm -D -u "$$file" | awk '$$NF ~ /^duckdb_/ {print "Unexpected DuckDB API import: " $$0; bad=1} END {exit bad}'; \
		fi; \
		test "$$actual" = duckvep_init_c_api || { printf 'Unexpected exports: %s\n' "$$actual"; exit 1; }
readme: release
	Rscript scripts/render-readme.R

.PHONY: check-rduckvep-bundle
check-rduckvep-bundle:
	scripts/check-rduckvep-bundle.sh

.PHONY: site site-check
site:
	Rscript scripts/build-site.R

site-check: site
	Rscript scripts/check-site-links.R

# AddressSanitizer + UBSan build of the extension in a separate directory, and the SQL suite
# run against it. The extension has no DuckDB link dependency (stable C API), so the runtime
# is preloaded into the unsanitized Python DuckDB host that the sqllogictest runner embeds.
# Any sanitizer report aborts the host process, which fails the target. Leak detection is off
# because the host process itself is not leak-clean.
ASAN_SANITIZE=-fsanitize=address,undefined -fno-omit-frame-pointer -fno-sanitize-recover=all
ASAN_LIB=$(shell gcc -print-file-name=libasan.so)
ASAN_CXXLIB=$(shell gcc -print-file-name=libstdc++.so.6)
build_asan: check_configure
	cmake $(CMAKE_BUILD_FLAGS) -DCMAKE_BUILD_TYPE=RelWithDebInfo \
		-DCMAKE_C_FLAGS="$(ASAN_SANITIZE) -O1 -g" -DCMAKE_CXX_FLAGS="$(ASAN_SANITIZE) -O1 -g" \
		-DCMAKE_SHARED_LINKER_FLAGS="$(ASAN_SANITIZE)" -S $(PROJ_DIR) -B cmake_build/asan
	cmake --build cmake_build/asan --config RelWithDebInfo
	mkdir -p build/asan
	$(PYTHON_VENV_BIN) extension-ci-tools/scripts/append_extension_metadata.py \
		-l cmake_build/asan/$(EXTENSION_LIB_FILENAME) -o build/asan/$(EXTENSION_FILENAME) \
		-n $(EXTENSION_NAME) -dv $(TARGET_DUCKDB_VERSION) \
		-evf configure/extension_version.txt -pf configure/platform.txt
test_release_asan: build_asan
	@test -f "$(ASAN_LIB)" || { echo "libasan.so not found (install gcc's libasan)" >&2; exit 1; }
	LD_PRELOAD="$(ASAN_LIB) $(ASAN_CXXLIB)" ASAN_OPTIONS=detect_leaks=0:abort_on_error=1:halt_on_error=1 \
		UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1 \
		$(TEST_RUNNER) --test-dir test/sql --external-extension build/asan/$(EXTENSION_FILENAME)
