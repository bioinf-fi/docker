#!/bin/bash
# Display installed genome assembly tools

check_tool() {
    local name="$1"
    local cmd="$2"
    printf "%s ...checking" "$name"
    local version=$(eval "$cmd" 2>&1)
    printf "\r\033[K%s %s\n" "$name" "$version"
}

echo "=========================================="
echo "Genome Assembly Tools"
echo "=========================================="
printf "%-12s %s\n" "Tool" "Version"
echo "------------------------------------------"

check_tool "verkko:" "verkko --version 2>&1 | head -n 1 || echo 'installed'"
check_tool "hifiasm:" "hifiasm --version 2>&1"
check_tool "seqtk:" "seqtk 2>&1 | grep -i version | head -n 1 | sed 's/Version: //'"
check_tool "samtools:" "samtools --version 2>&1 | head -n 1 | cut -d' ' -f2"
check_tool "minimap2:" "minimap2 --version 2>&1"
check_tool "mashmap:" "mashmap --version 2>&1"
check_tool "bioawk:" "bioawk --version 2>&1 | head -n 1 | awk '{print \$3}'"
check_tool "quast:" "quast --version 2>&1 | grep -i quast | sed 's/QUAST //'"
check_tool "nanoplot:" "python -c 'import importlib.metadata; print(importlib.metadata.version(\"NanoPlot\"))' 2>/dev/null || echo 'installed'"

echo ""
echo "GUI Tools (use 'make gui' to run):"
echo "------------------------------------------"
check_tool "IGV:" "igv --version 2>&1 | head -n 1 || echo 'installed'"
check_tool "Bandage:" "Bandage --version 2>&1 | grep -i version | head -n 1 || echo 'installed'"

echo "=========================================="
echo "Working directory: /data"
echo "=========================================="
