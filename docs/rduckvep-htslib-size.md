# Indexed FASTA extension size

These Linux amd64 measurements compare two build configurations, not parent and child revisions:

| Build | Bytes | Provenance |
| --- | ---: | --- |
| Full HTSlib, isolated design study | 1,397,206 | HTSlib design memo, revision `2d613d92` |
| zlib-only HTSlib with section garbage collection | 609,614 | `make release` at commit `d3babd4`; `stat -c %s build/release/duckvep.duckdb_extension` |

Both sizes include DuckDB extension metadata. The 609,614-byte Linux shared library dynamically links zlib, libm and libc; it does not link curl, crypto, bzip2, lzma or libdeflate. At `d3babd4`, `make test_release` passed all eight SQL tests and the SQL lambda guard. Windows and macOS dead-strip measurements require platform-specific builds.
