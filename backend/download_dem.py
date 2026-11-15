"""Download the NASA/USGS MOLA DEM (about 2 GB, resumable) into backend/data/."""
import argparse
import os
import sys
import urllib.request

URL = "https://planetarymaps.usgs.gov/mosaic/Mars_MGS_MOLA_DEM_mosaic_global_463m.tif"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "Mars_MGS_MOLA_DEM_mosaic_global_463m.tif")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--url", default=URL)
    ap.add_argument("--out", default=OUT)
    args = ap.parse_args()
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    part = args.out + ".part"
    have = os.path.getsize(part) if os.path.exists(part) else 0
    headers = {
    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 Chrome/120.0 Safari/537.36",
    "Accept": "*/*"
    }

    if have:
        headers["Range"] = f"bytes={have}-"

    req = urllib.request.Request(args.url, headers=headers)
    with urllib.request.urlopen(req) as resp:
        total = have + int(resp.headers.get("Content-Length", 0))
        if have and resp.status != 206:  # server ignored Range: start over
            have = 0
        with open(part, "ab" if have else "wb") as fh:
            done = have
            while chunk := resp.read(1 << 20):
                fh.write(chunk)
                done += len(chunk)
                sys.stdout.write(f"\r{done/1e6:8.0f} / {total/1e6:.0f} MB")
                sys.stdout.flush()
    os.replace(part, args.out)
    print(f"\nSaved {args.out}")


if __name__ == "__main__":
    main()
