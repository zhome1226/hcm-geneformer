"""Download author-supplied SCP1303 files from the official HCA S3 mirror.

Only transport and integrity checks; does not load expression matrices.
Input is a saved response from HCA's index/files API for project
85c0d6fa-f117-4d76-b01a-5d5e8f5f9188. Standard library + curl required.
"""
import argparse
import concurrent.futures
import csv
import gzip
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
import shutil
import subprocess
import time

SELECTED = {
    "DCM_HCM_Expression_Matrix_barcodes_V1.tsv",
    "DCM_HCM_Expression_Matrix_genes_V1.tsv",
    "DCM_HCM_MetaData_V1.txt",
    "DCM_HCM_UMAP_V1.txt",
    "human_dcm_hcm_scportal_03.17.2022.h5ad.gz",
    "human_dcm_hcm_scportal_03.17.2022_README.txt",
    "ChaffinHypertrophicCardiomyopathy_metadata_01-12-2023.xlsx",
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", required=True, type=Path)
    parser.add_argument("--base", required=True, type=Path)
    parser.add_argument("--proxy")
    args = parser.parse_args()
    base = args.base.resolve()
    base.mkdir(parents=True, exist_ok=True)
    logs = base / "_catalog"
    logs.mkdir(exist_ok=True)
    catalog = json.loads(args.catalog.read_text())
    files = [f for hit in catalog["hits"] for f in hit["files"] if f["name"] in SELECTED]
    if len(files) != len(SELECTED):
        raise ValueError("Official catalog does not contain every selected filename")
    needed = sum(max(0, f["size"] - (base / f["name"]).stat().st_size)
                 if (base / f["name"]).exists() else f["size"] for f in files)
    if shutil.disk_usage(base).free < needed + 10 * 1024**3:
        raise RuntimeError("Insufficient space: require payload plus 10 GiB reserve")

    def download(f):
        name = f["name"]
        target = base / name
        url = "https://s3.amazonaws.com/humancellatlas/file/" + f["sha256"] + ".sha256"
        print("START " + name, flush=True)
        if not target.exists() or target.stat().st_size != f["size"]:
            command = ["curl", "--fail", "--location", "--continue-at", "-",
                       "--retry", "0", "--connect-timeout", "20",
                       "--speed-limit", "1024", "--speed-time", "180",
                       "--output", str(target)]
            if args.proxy:
                command.extend(["--proxy", args.proxy])
            command.append(url)
            with (logs / (name + ".download.log")).open("ab") as log:
                # Launch a fresh curl on every retry so -C - recomputes the
                # actual file length. Curl's internal retries can reset to the
                # offset captured at the beginning of its first attempt.
                for attempt in range(9):
                    result = subprocess.run(command, stdout=log, stderr=log)
                    if result.returncode == 0:
                        break
                    if attempt == 8:
                        result.check_returncode()
                    print(f"RESUME {name} at {target.stat().st_size if target.exists() else 0} bytes", flush=True)
                    time.sleep(min(2 ** attempt, 30))
        if target.stat().st_size != f["size"]:
            raise ValueError("Size mismatch: " + name)
        sha256, md5 = hashlib.sha256(), hashlib.md5()
        with target.open("rb") as source:
            for block in iter(lambda: source.read(8 * 1024**2), b""):
                sha256.update(block)
                md5.update(block)
        if sha256.hexdigest() != f["sha256"]:
            raise ValueError("Official SHA-256 mismatch: " + name)
        if name.endswith(".gz"):
            with gzip.open(target, "rb") as source:
                while source.read(8 * 1024**2):
                    pass
        row = {"dataset": "SCP1303", "file": name, "relative_path": name,
               "bytes": target.stat().st_size, "md5": md5.hexdigest(),
               "sha256": sha256.hexdigest(), "source_url": url,
               "verification": "official HCA SHA-256 and size match",
               "disease_groups": "15 HCM;11 DCM;16 NF",
               "tissue": "human left ventricle", "purpose": "planned training source",
               "download_date_utc": datetime.now(timezone.utc).isoformat()}
        (logs / (name + ".verified.json")).write_text(json.dumps(row, indent=2) + "\n")
        print("VERIFIED " + name, flush=True)
        return row

    rows = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        for row in pool.map(download, files):
            rows.append(row)
    rows.sort(key=lambda row: row["file"])
    with (base / "MANIFEST.tsv").open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    (base / "MD5SUMS").write_text("".join(r["md5"] + "  " + r["file"] + "\n" for r in rows))
    (base / "SHA256SUMS").write_text("".join(r["sha256"] + "  " + r["file"] + "\n" for r in rows))
    summary = {"status": "complete", "files": len(rows), "bytes": sum(r["bytes"] for r in rows),
               "completed_utc": datetime.now(timezone.utc).isoformat(), "files_verified": rows}
    (logs / "download_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "files_verified"}), flush=True)


if __name__ == "__main__":
    main()
