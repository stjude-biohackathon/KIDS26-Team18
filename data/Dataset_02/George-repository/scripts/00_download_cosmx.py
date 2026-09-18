#!/usr/bin/env python
from pathlib import Path
from urllib.request import urlopen, Request
import shutil, sys, time
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from src.cosmx_io import BASE_URL, FILES, validate_gzip_csv
RAW=Path(__file__).resolve().parents[1]/"data/raw"; RAW.mkdir(parents=True,exist_ok=True)
for kind,name in FILES.items():
    dst=RAW/name
    if dst.exists():
        try: validate_gzip_csv(dst); print(f"[skip] valid: {dst.name}"); continue
        except Exception: print(f"[redo] invalid existing file: {dst.name}")
    url=BASE_URL+name; tmp=dst.with_suffix(dst.suffix+".part")
    print(f"[download] {kind}: {url}")
    for attempt in range(1,4):
        try:
            req=Request(url,headers={"User-Agent":"Mozilla/5.0"})
            with urlopen(req,timeout=120) as r, open(tmp,"wb") as w: shutil.copyfileobj(r,w,length=8*1024*1024)
            tmp.replace(dst); validate_gzip_csv(dst); break
        except Exception as e:
            print(f"  attempt {attempt}/3 failed: {e}")
            if tmp.exists(): tmp.unlink()
            if attempt==3: raise
            time.sleep(5*attempt)
print("All five GEO CosMx flat files are present and readable.")
