import glob, sys

wheels = glob.glob("/wheelhouse/*.whl")
need = {"pysam","moddotplot","whatshap","pybedtools","pybigwig","ndindex","methylartist"}
have = {w.split("/")[-1].split("-")[0].lower() for w in wheels}
missing = sorted(need - have)
if missing:
    print("Missing wheels:", missing)
    sys.exit(1)
print("Wheelhouse OK.")

