#!/bin/bash
pip install -e .

echo "Building CaloX utilities..."
root -l -q -b -e '.L utils/functions.cc+' && echo "✓ Build complete"
