# batch_run_curv_methods.py
# Direct run, no if __name__ == "__main__"
from pathlib import Path
import re
import shutil
import subprocess
from datetime import datetime

# ===================== USER CONFIG =====================
TEMPLATE_IN = Path(r"./cur.inp")  
EXE = Path(r"/home/mingrgao/sisc_lab/fused_backup/curv").resolve()

OUT_ROOT    = Path(r"./batch_runs")        

METHODS = ["CV", "RDF_raw", "RDF_smooth", "PARABOLOID", "HF", "VFM", "PC", "PV"]


OVERRIDE_N_LIST = None  # e.g. [20, 40, 80, 120, 160]
# =======================================================


def read_text(p: Path) -> str:
    return p.read_text(encoding="utf-8", errors="ignore")

def write_text(p: Path, s: str):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(s, encoding="utf-8")

def patch_curv_method(text: str, method: str) -> str:
    pat = r"(?m)^\s*curv\.method\s*=\s*\S+\s*$"
    repl = f"curv.method = {method}"
    if re.search(pat, text):
        return re.sub(pat, repl, text, count=1)
    return text.rstrip() + "\n\n" + repl + "\n"

def patch_n_list(text: str, n_list: list[int]) -> str:
    pat = r"(?m)^\s*amr\.n_list\s*=\s*.+\s*$"
    repl = "amr.n_list = " + " ".join(str(x) for x in n_list)
    if re.search(pat, text):
        return re.sub(pat, repl, text, count=1)
    return text.rstrip() + "\n\n" + repl + "\n"

def parse_n_list(text: str) -> list[int]:
    m = re.search(r"(?m)^\s*amr\.n_list\s*=\s*(.+)\s*$", text)
    if not m:
        return []
    out = []
    for tok in m.group(1).split():
        tok = tok.strip()
        if tok.isdigit():
            out.append(int(tok))
    return out

def run_try_arg(exe: Path, input_file: Path, workdir: Path, log_file: Path) -> int:
   
    cmd = [str(exe), str(input_file.name)]
    with log_file.open("wb") as fout:
        p = subprocess.run(cmd, stdout=fout, stderr=subprocess.STDOUT, cwd=str(workdir))
    return p.returncode

def run_try_stdin(exe: Path, input_file: Path, workdir: Path, log_file: Path) -> int:
  
    cmd = [str(exe)]
    with input_file.open("rb") as fin, log_file.open("wb") as fout:
        p = subprocess.run(cmd, stdin=fin, stdout=fout, stderr=subprocess.STDOUT, cwd=str(workdir))
    return p.returncode

def looks_like_input_open_error(log_text: str) -> bool:
    t = log_text.lower()
    hints = [
        "usage", "argument", "no such file", "cannot open", "failed to open",
        "input", "filename", "expects", "missing"
    ]
    return any(h in t for h in hints)

def run_auto(exe: Path, input_file: Path, workdir: Path) -> tuple[int, str]:
    log1 = workdir / "run_arg.log"
    rc1 = run_try_arg(exe, input_file, workdir, log1)

    if rc1 == 0:
       
        (workdir / "run.log").write_bytes(log1.read_bytes())
        log1.unlink(missing_ok=True)
        return 0, "arg"

    txt1 = ""
    try:
        txt1 = log1.read_text(encoding="utf-8", errors="ignore")
    except Exception:
        pass

    if looks_like_input_open_error(txt1):
        log2 = workdir / "run_stdin.log"
        rc2 = run_try_stdin(exe, input_file, workdir, log2)

        
        if rc2 == 0:
            (workdir / "run.log").write_bytes(log2.read_bytes())
            log1.unlink(missing_ok=True)
            log2.unlink(missing_ok=True)
            return 0, "stdin"
        else:
            
            (workdir / "run.log").write_bytes(log1.read_bytes())
            return rc1, "arg_failed_stdin_failed"

    (workdir / "run.log").write_bytes(log1.read_bytes())
    log1.unlink(missing_ok=True)
    return rc1, "arg"


# ---------------- main flow ----------------
if not TEMPLATE_IN.exists():
    raise FileNotFoundError(f"TEMPLATE_IN not found: {TEMPLATE_IN}")
if not EXE.exists():
    raise FileNotFoundError(f"EXE not found: {EXE}")

tmpl = read_text(TEMPLATE_IN)

n_list = OVERRIDE_N_LIST if OVERRIDE_N_LIST is not None else parse_n_list(tmpl)
if not n_list:
    n_list = [None]  

stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
root = OUT_ROOT / f"run_{stamp}"
root.mkdir(parents=True, exist_ok=True)
shutil.copy2(TEMPLATE_IN, root / "inputs_template.in")

summary = []
summary.append(f"root: {root}")
summary.append(f"exe : {EXE.resolve()}")
summary.append(f"methods: {', '.join(METHODS)}")
summary.append(f"n_list: {n_list}")

for method in METHODS:
    for n in n_list:
        ntag = f"n{n:04d}" if isinstance(n, int) else "nNA"
        outdir = root / method / ntag
        outdir.mkdir(parents=True, exist_ok=True)

        text = patch_curv_method(tmpl, method)
        if isinstance(n, int):
            text = patch_n_list(text, [n])

        input_path = outdir / "inputs.in"
        write_text(input_path, text)

        rc, mode = run_auto(EXE, input_path, outdir)
        summary.append(f"{method:10s} {ntag:6s} rc={rc} mode={mode}  log={outdir.relative_to(root) / 'run.log'}")

write_text(root / "summary.txt", "\n".join(summary) + "\n")

print("\n".join(summary))
print(f"\nDONE. Summary: {root / 'summary.txt'}")
