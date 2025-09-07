import importlib, sys

mods = (
    "pysam",
    "whatshap",
    "moddotplot",
    "pybedtools",
    "pyBigWig",
    "ndindex",
    "modbamtools",
)
bad = []
for m in mods:
    try:
        importlib.import_module(m)
    except Exception as e:
        bad.append((m, e.__class__.__name__, str(e)))
if bad:
    print("Import failures:", bad)
    sys.exit(1)
print("Core Python modules import OK.")
