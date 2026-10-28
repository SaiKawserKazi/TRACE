"""Download every input listed in DATA_MANIFEST.csv into data/raw/ and check its SHA-256 (cross-platform; resumes partial files).

    python src/download.py            # download missing files, verify all
    python src/download.py --check    # verify only

Some sources are live lists without a version (OncoKB, Enrichr Hallmark, HAGR, Alliance, BAGEL CEGv2); for these a checksum
mismatch is reported as a warning, not an error, because the provider may have updated the file since our access date."""
import csv, hashlib, sys, urllib.request
from config import ROOT, RAW

LIVE = ("oncokb", "Enrichr", "genomics.senescence.info", "alliancegenome", "hart-lab")


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()


def fetch(url, dest):
    part = dest.with_name(dest.name + ".part"); start = part.stat().st_size if part.exists() else 0
    req = urllib.request.Request(url, headers={"User-Agent": "TRACE-download/1.0", **({"Range": f"bytes={start}-"} if start else {})})
    with urllib.request.urlopen(req, timeout=120) as r, open(part, "ab" if start and r.status == 206 else "wb") as f:
        while True:
            b = r.read(1 << 20)
            if not b: break
            f.write(b)
    part.replace(dest)


def main():
    check_only = "--check" in sys.argv; bad = 0
    for row in csv.DictReader(open(ROOT / "DATA_MANIFEST.csv", encoding="utf-8-sig")):
        dest = RAW / row["local_file"]
        if not dest.exists():
            if check_only: print("MISSING", dest.name); bad += 1; continue
            print("downloading", dest.name, "from", row["url"], flush=True); fetch(row["url"], dest)
        ok = sha256(dest) == row["sha256"]
        live = any(k in row["url"] for k in LIVE)
        print(("ok      " if ok else ("CHANGED " if live else "MISMATCH")), dest.name)
        bad += (not ok) and not live
    if bad: sys.exit(f"{bad} file(s) missing or with a wrong checksum")


if __name__ == "__main__":
    main()
