# Genomics Toolkit helpers

tools() {
  echo "Tools ready:"
  command -v samtools >/dev/null 2>&1 \
    && samtools --version 2>/dev/null | head -n1 \
    || echo "samtools: not found (or missing libhts)"
  # minimap2 (print name + version)
  if command -v minimap2 >/dev/null 2>&1; then
    v="$(minimap2 --version 2>/dev/null || true)"
    [ -n "$v" ] && echo "minimap2 $v" || echo "minimap2: installed"
  else
    echo "minimap2: not found"
  fi
  # seqtk (print name + version)
    if command -v seqtk >/dev/null 2>&1; then
    ver="$(seqtk 2>&1 | awk '/^Version:/{v=$2; if ($3) v=v"-"$3; print v; exit}')"
    if [ -n "$ver" ]; then
    echo "seqtk $ver"
    else
    # fallback: first non-empty line of help
    line="$(seqtk 2>&1 | awk 'NF{print;exit}')"
    [ -n "$line" ] && echo "seqtk $line" || echo "seqtk: installed"
    fi
    else
    echo "seqtk: not found"
    fi
  command -v bioawk >/dev/null 2>&1 \
    && echo "bioawk: installed" \
    || echo "bioawk: not found"

  # Import-check for core Python libs (exclude methylartist)
  /opt/venv/bin/python - <<'PY'
import importlib
mods=("pysam","whatshap","moddotplot","pybedtools","pyBigWig","ndindex")
for m in mods:
    try:
        mod=importlib.import_module(m)
        print(f"{m} {getattr(mod,'__version__','installed')}")
    except Exception as e:
        print(f"{m}: import failed ({e.__class__.__name__}: {e})")
PY

  # Methylartist via CLI
  if command -v methylartist >/dev/null 2>&1; then
    # some versions print help on --version; fall back to a simple presence message
    methylartist --version 2>/dev/null || echo "methylartist: installed (CLI)"
  else
    echo "methylartist: not found"
  fi
