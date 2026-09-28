# Indexed FASTA build size

| Linux amd64 release artifact | Bytes | Source |
| --- | ---: | --- |
| Before (full HTSlib) | 1,397,206 | HTSlib design memo, revision `2d613d92` |
| After (zlib-only, section GC) | 609,614 | `make release`, this checkout, `stat -c %s build/release/duckvep.duckdb_extension` |

The before measurement was made in the isolated design study, not by building the parent of this commit; toolchain and revision differ. Both numbers include DuckDB extension metadata. The resulting Linux shared library dynamically links zlib, libm and libc, not curl, crypto, bzip2, lzma or libdeflate. `make test_release` passed all eight SQL tests and the SQL lambda guard. Windows and macOS dead-strip results require their own builds.
