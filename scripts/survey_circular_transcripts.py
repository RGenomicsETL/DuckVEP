#!/usr/bin/env python3
"""Survey public Ensembl core databases for origin-crossing transcripts.

For every core database under an Ensembl MySQL dump directory this reads the small
attrib_type and seq_region_attrib tables to find regions with the circular_seq
attribute, then, only for databases that have one, streams transcript.txt.gz and
counts transcripts on those regions with seq_region_start > seq_region_end. Nothing
is stored beyond the printed table; tables are streamed and decoded on the fly.

    scripts/survey_circular_transcripts.py https://ftp.ensembl.org/pub/release-116/mysql/ \
        https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-63/plants/mysql/ ...
"""
import gzip
import re
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor


def open_table(base, database, table, attempts=6):
    """Open one gzip table, retrying: the FTP servers refuse bursts of connections."""
    url = f"{base}{database}/{table}.txt.gz"
    for attempt in range(attempts):
        try:
            request = urllib.request.Request(url)
            return gzip.open(urllib.request.urlopen(request, timeout=120), "rt", encoding="utf-8",
                             errors="replace")
        except OSError:
            if attempt == attempts - 1:
                raise
            time.sleep(2 ** attempt)


def databases(base):
    with urllib.request.urlopen(base, timeout=120) as handle:
        html = handle.read().decode()
    return sorted(set(re.findall(r'href="([a-z0-9_]+_core_[0-9_]+)/"', html)))


def survey(item):
    base, database = item
    try:
        code_id = None
        with open_table(base, database, "attrib_type") as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) > 1 and fields[1] == "circular_seq":
                    code_id = fields[0]
                    break
        if code_id is None:
            return database, None, 0, 0, 0
        regions = set()
        with open_table(base, database, "seq_region_attrib") as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) > 2 and fields[1] == code_id and fields[2] == "1":
                    regions.add(fields[0])
        if not regions:
            return database, None, 0, 0, 0
        transcripts = wrapped = 0
        with open_table(base, database, "transcript") as handle:
            for line in handle:
                fields = line.split("\t", 7)
                if fields[3] in regions:
                    transcripts += 1
                    if int(fields[4]) > int(fields[5]):
                        wrapped += 1
        return database, len(regions), transcripts, wrapped, 0
    except Exception as error:  # a missing table is reported, not fatal
        return database, "error: %s" % error, 0, 0, 0


def main():
    items = [(base, database) for base in sys.argv[1:] for database in databases(base)]
    print("core databases: %d" % len(items), flush=True)
    with ThreadPoolExecutor(max_workers=3) as pool:
        rows = list(pool.map(survey, items))
    circular = [r for r in rows if isinstance(r[1], int)]
    errors = [r for r in rows if isinstance(r[1], str)]
    print("databases with circular_seq=1 regions: %d" % len(circular))
    print("databases that failed: %d" % len(errors))
    print("database\tcircular_regions\ttranscripts_on_them\tinverted_transcripts")
    for database, regions, transcripts, wrapped, _ in circular:
        print("%s\t%d\t%d\t%d" % (database, regions, transcripts, wrapped))
    for database, message, *_ in errors:
        print("ERROR\t%s\t%s" % (database, message))
    print("total inverted transcripts: %d" % sum(r[3] for r in circular))


if __name__ == "__main__":
    main()
