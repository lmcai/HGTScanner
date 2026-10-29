#!/usr/bin/env python3
"""
hgtscanner_config.py - set up and update the sequence databases for HGTScanner.

The large FASTA databases are hosted on Zenodo (gzipped). This script downloads
them, verifies the MD5 checksums, and unpacks them into <HGTScanner>/database,
where HGTscanner.py finds them automatically.

Place this script in the top-level HGTScanner folder (next to HGTscanner.py).
Only the Python standard library is needed (Python >= 3.6).

Usage:
    python hgtscanner_config.py            # first-time setup: download missing databases
    python hgtscanner_config.py --update   # replace the databases with the newest release

Database file names carry a release date (e.g. Viridiplantae_mt_aug2025.genome.fas).
HGTscanner.py always picks the newest release present in the database folder, so
no source code edits are needed after an update.
"""

import argparse
import gzip
import hashlib
import json
import re
import shutil
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

ZENODO_RECORD = "23267980"          # any version of the Zenodo record works; the newest is resolved automatically
API_BASE = "https://zenodo.org/api/records"
DB_DIR = Path(__file__).resolve().parent / "database"

CHUNK = 1024 * 1024
RETRIES = 5
USER_AGENT = "HGTScanner-config/1.0"

# <prefix>_<monYYYY>.<kind>.fas  ->  e.g. Viridiplantae_mt_aug2025.genome.fas
DATED = re.compile(r"^(?P<prefix>.+?)_(?P<date>[a-z]{3}\d{4})\.(?P<kind>[^.]+)\.(?:fas|fasta|fa)$", re.I)


# ----------------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------------
def human(n):
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024 or unit == "TB":
            return f"{int(n)} B" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024


def http_open(url, headers=None, timeout=60):
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, **(headers or {})})
    return urllib.request.urlopen(req, timeout=timeout)


def fetch_json(url):
    last = None
    for attempt in range(1, RETRIES + 1):
        try:
            with http_open(url) as r:
                return json.load(r)
        except urllib.error.HTTPError as e:
            if e.code == 404:
                raise
            last = e
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as e:
            last = e
        time.sleep(2 * attempt)
    raise last


def md5sum(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(CHUNK), b""):
            h.update(block)
    return h.hexdigest()


def key_of(fasta_name):
    """Identity of a database independent of its release date."""
    m = DATED.match(fasta_name)
    return (m.group("prefix"), m.group("kind")) if m else (fasta_name, None)


def get_latest_record(record):
    """Return the metadata of the newest version of the Zenodo record."""
    for url in (f"{API_BASE}/{record}/versions/latest", f"{API_BASE}/{record}"):
        try:
            return fetch_json(url)
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError):
            continue
    sys.exit(f"ERROR: could not reach Zenodo record {record}. Check your internet connection.")


def list_remote_files(meta):
    files = []
    for f in meta.get("files", []):
        name = f.get("key") or f.get("filename")
        links = f.get("links", {})
        url = links.get("self") or links.get("download") or \
            f"https://zenodo.org/records/{meta.get('id')}/files/{name}?download=1"
        checksum = f.get("checksum", "")
        md5 = checksum.split(":", 1)[1] if checksum.startswith("md5:") else (checksum or None)
        files.append({"name": name, "url": url, "size": f.get("size", 0), "md5": md5})
    if not files:
        sys.exit("ERROR: no files found in the Zenodo record.")
    return files


# ----------------------------------------------------------------------------
# download (resumable) and decompress
# ----------------------------------------------------------------------------
def download(url, dest, size, md5):
    part = dest.with_name(dest.name + ".part")
    for attempt in range(1, RETRIES + 1):
        try:
            have = part.stat().st_size if part.exists() else 0
            if size and have > size:
                part.unlink()
                have = 0
            headers = {"Range": f"bytes={have}-"} if have else {}
            start = time.time()
            with http_open(url, headers=headers) as r:
                if have and r.status != 206:      # server ignored Range -> start over
                    have = 0
                total = size or (int(r.headers.get("Content-Length", 0)) + have)
                done = have
                with open(part, "ab" if have else "wb") as out:
                    while True:
                        block = r.read(CHUNK)
                        if not block:
                            break
                        out.write(block)
                        done += len(block)
                        speed = (done - have) / max(time.time() - start, 1e-6)
                        pct = f"{100 * done / total:5.1f}%  " if total else ""
                        sys.stdout.write(f"\r    {pct}{human(done)} / {human(total)}  ({human(speed)}/s)    ")
                        sys.stdout.flush()
            print()
            if size and part.stat().st_size != size:
                raise IOError(f"size mismatch ({part.stat().st_size} != {size})")
            if md5:
                sys.stdout.write("    verifying MD5 ... ")
                sys.stdout.flush()
                if md5sum(part) != md5:
                    part.unlink()
                    raise IOError("MD5 checksum mismatch")
                print("OK")
            part.replace(dest)
            return
        except (urllib.error.URLError, TimeoutError, IOError, ConnectionError) as e:
            print(f"\n    download problem: {e}")
            if attempt == RETRIES:
                sys.exit(f"ERROR: failed to download {dest.name} after {RETRIES} attempts. Re-run to resume.")
            wait = 3 * attempt
            print(f"    retrying in {wait}s (attempt {attempt + 1}/{RETRIES}) ...")
            time.sleep(wait)


def gunzip(src, dest):
    tmp = dest.with_name(dest.name + ".tmp")
    start = time.time()
    written = 0
    with gzip.open(src, "rb") as fin, open(tmp, "wb") as fout:
        while True:
            block = fin.read(CHUNK)
            if not block:
                break
            fout.write(block)
            written += len(block)
            speed = written / max(time.time() - start, 1e-6)
            sys.stdout.write(f"\r    unpacking: {human(written)}  ({human(speed)}/s)    ")
            sys.stdout.flush()
    print()
    tmp.replace(dest)


# ----------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Download and update the HGTScanner sequence databases from Zenodo.")
    ap.add_argument("--update", action="store_true",
                    help="download the newest database release and replace the old FASTA files")
    ap.add_argument("--record", default=ZENODO_RECORD,
                    help=f"Zenodo record ID (default: {ZENODO_RECORD}; the newest version is used)")
    args = ap.parse_args()

    DB_DIR.mkdir(parents=True, exist_ok=True)
    print(f"HGTScanner database {'update' if args.update else 'setup'}\n  Target folder: {DB_DIR}\n")

    meta = get_latest_record(args.record)
    remote = list_remote_files(meta)
    print(f"Newest Zenodo release: record {meta.get('id')} "
          f"(published {meta.get('metadata', {}).get('publication_date', 'n/a')}), "
          f"{len(remote)} file(s), {human(sum(f['size'] for f in remote))} compressed.")
    print(f"Free disk space: {human(shutil.disk_usage(DB_DIR).free)} "
          f"(unpacked FASTA files need several times the download size)\n")

    local = [p for p in DB_DIR.iterdir() if p.is_file() and DATED.match(p.name)]
    newer_available = False

    for i, f in enumerate(remote, 1):
        is_gz = f["name"].endswith(".gz")
        fasta_name = f["name"][:-3] if is_gz else f["name"]
        target = DB_DIR / fasta_name
        gz_path = DB_DIR / f["name"]
        key = key_of(fasta_name)
        old = [p for p in local if key[1] and key_of(p.name) == key and p.name != fasta_name]

        print(f"[{i}/{len(remote)}] {fasta_name}  ({human(f['size'])})")

        if target.exists() and target.stat().st_size > 0:
            print("    up to date\n")
            continue
        if old and not args.update:
            print(f"    newer release of {', '.join(p.name for p in old)} is available "
                  f"(not downloaded; use --update)\n")
            newer_available = True
            continue

        if gz_path.exists() and f["md5"] and md5sum(gz_path) == f["md5"]:
            print("    archive already downloaded and verified")
        else:
            download(f["url"], gz_path, f["size"], f["md5"])
        if is_gz:
            gunzip(gz_path, target)
            gz_path.unlink()
        for p in old:                       # replace old release only after the new one is in place
            p.unlink()
            print(f"    removed old release {p.name}")
        print(f"    -> {target.name}\n")

    print("Database folder contents:")
    for p in sorted(DB_DIR.iterdir()):
        if p.is_file() and p.suffix in (".fas", ".fasta", ".fa"):
            print(f"  {p.name:50s} {human(p.stat().st_size)}")
    if newer_available:
        print("\nA newer database release exists. Run:  python hgtscanner_config.py --update")
    else:
        print("\nDatabase is ready. You can now run HGTscanner.py.")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.exit("\nInterrupted. Re-run the script to resume.")
