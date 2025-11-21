#!/bin/bash
# VCA Pipeline Launcher

set -e

CONFIG_FILE="pipeline_config.yml"

if [ $# -eq 1 ]; then
    CONFIG_FILE="$1"
fi

if [ ! -f "$CONFIG_FILE" ]; then
    echo "[ERROR] Config file not found: ${CONFIG_FILE}"
    exit 1
fi

if [ -z "$CONDA_DEFAULT_ENV" ]; then
    echo "[WARNING] Conda environment not activated"
    echo "Please run: conda activate vca_env"
    exit 1
fi

echo "Checking dependencies..."
for tool in bwa samtools bcftools prefetch fasterq-dump; do
    if ! command -v $tool &> /dev/null; then
        echo "[ERROR] Required tool not found: $tool"
        exit 1
    fi
done
echo "[OK] All dependencies found"

echo ""
echo "VCA Pipeline Configuration"
echo "Config: ${CONFIG_FILE}"
echo "Conda: ${CONDA_DEFAULT_ENV}"
echo "Working dir: $(pwd)"
echo ""

read -p "Start pipeline? (y/n) " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Pipeline cancelled"
    exit 0
fi

echo "Starting VCA pipeline..."
python run_vca_pipeline.py "$CONFIG_FILE"

EXIT_STATUS=$?

if [ $EXIT_STATUS -eq 0 ]; then
    echo ""
    echo "Pipeline completed successfully!"
    echo "Results: ../data_cluster/results/"
else
    echo ""
    echo "Pipeline failed with exit code: ${EXIT_STATUS}"
fi

exit $EXIT_STATUS
