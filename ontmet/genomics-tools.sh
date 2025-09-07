# Genomics Toolkit helpers

tools() {
  echo "Tools ready:"
  command -v samtools >/dev/null 2>&1 \
    && samtools --version 2>/dev/null | head -n1 \
    || echo "samtools: not found (or missing libhts)"
  command -v bcftools >/dev/null 2>&1 \
    && bcftools --version 2>/dev/null | head -n1 \
    || echo "bcftools: not found (or missing libhts)"
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

  # Additional tools with version info
  if command -v modkit >/dev/null 2>&1; then
    modkit --version | head -1
  else
    echo "modkit: not found"
  fi

  if command -v igv >/dev/null 2>&1; then
    # IGV version is in the JAR manifest - extract from directory name or use generic
    igv_dir="$(dirname "$(readlink -f "$(command -v igv)")")"
    igv_ver="$(basename "$igv_dir" | grep -o '[0-9]\+\.[0-9]\+\.[0-9]\+' || echo "2.18.2")"
    echo "IGV ${igv_ver}"
  else
    echo "IGV: not found"
  fi

  if command -v pomfret >/dev/null 2>&1; then
    # pomfret doesn't have --version, but we can get basic info
    pomfret --help 2>&1 | head -1 | sed 's/^/pomfret: /' || echo "pomfret: installed"
  else
    echo "pomfret: not found"
  fi

  if command -v NanoPlot >/dev/null 2>&1; then
    NanoPlot --version 2>/dev/null | head -1 || echo "NanoPlot: installed"
  else
    echo "NanoPlot: not found"
  fi

  # Check for modbamtools (pip-installed in main venv)
  /opt/venv/bin/python -c "
try:
    import modbamtools
    # Try multiple ways to get version
    version = getattr(modbamtools, '__version__', None)
    if not version:
        # Try from package metadata
        try:
            import importlib.metadata
            version = importlib.metadata.version('modbamtools')
        except:
            # Try from pkg_resources (older method)
            try:
                import pkg_resources
                version = pkg_resources.get_distribution('modbamtools').version
            except:
                version = 'installed'
    print(f'modbamtools {version}')
except ImportError:
    print('modbamtools: not found')
except Exception as e:
    print(f'modbamtools: error ({e})')
" 2>/dev/null
}
