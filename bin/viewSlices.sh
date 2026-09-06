#!/bin/bash
# VIEWSLICES: Launch Fast Tractography Slice Viewer (Local Computer)
#
# This script launches the Python viewer on your local computer.
# NO MATLAB REQUIRED - only reads pre-generated slice images.
#
# Usage:
#   ./viewSlices.sh /path/to/cache
#   ./viewSlices.sh                    (will prompt for directory)
#
# Prerequisites:
#   • Python 3.7+
#   • pip install pillow numpy

VIEWER_SCRIPT="FastTractographyViewer.py"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "╔═══════════════════════════════════════════════════════════╗"
echo "║   Fast Tractography Slice Viewer (Local Computer)        ║"
echo "╚═══════════════════════════════════════════════════════════╝"
echo ""

# Check if Python is available
if ! command -v python3 &> /dev/null; then
    if ! command -v python &> /dev/null; then
        echo "✗ ERROR: Python not found"
        echo "  Please install Python 3.7+ from python.org"
        exit 1
    fi
    PYTHON_CMD="python"
else
    PYTHON_CMD="python3"
fi

echo "✓ Python: $($PYTHON_CMD --version)"

# Check if viewer script exists
# The viewer lives in scripts/, not beside this launcher in bin/. Resolving it
# against SCRIPT_DIR made this script always fail with "not found".
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
VIEWER_PATH="$REPO_ROOT/scripts/$VIEWER_SCRIPT"
if [ ! -f "$VIEWER_PATH" ]; then
    echo "✗ ERROR: $VIEWER_SCRIPT not found"
    echo "  Expected location: $VIEWER_PATH"
    exit 1
fi

echo "✓ Viewer: $VIEWER_SCRIPT"

# Check dependencies
echo ""
echo "Checking dependencies..."

check_package() {
    if $PYTHON_CMD -c "import $1" 2>/dev/null; then
        echo "  ✓ $1"
        return 0
    else
        echo "  ✗ $1 (missing)"
        return 1
    fi
}

MISSING_DEPS=0
check_package "tkinter" || MISSING_DEPS=1
check_package "PIL" || MISSING_DEPS=1
check_package "numpy" || MISSING_DEPS=1

if [ $MISSING_DEPS -eq 1 ]; then
    echo ""
    echo "Installing missing dependencies..."
    $PYTHON_CMD -m pip install pillow numpy
    echo ""
fi

# Get cache directory
if [ -z "$1" ]; then
    echo ""
    echo "Enter path to slice cache directory:"
    echo "(The directory transferred from the server)"
    read -e -p "> " CACHE_DIR
else
    CACHE_DIR="$1"
fi

# Validate cache directory
if [ ! -d "$CACHE_DIR" ]; then
    echo ""
    echo "✗ ERROR: Directory not found: $CACHE_DIR"
    exit 1
fi

if [ ! -d "$CACHE_DIR/datasets" ]; then
    echo ""
    echo "✗ ERROR: Invalid cache directory (missing datasets/)"
    echo "  Make sure you're pointing to the root cache directory"
    echo "  transferred from the server."
    exit 1
fi

echo ""
echo "✓ Cache directory: $CACHE_DIR"
echo ""
echo "═══════════════════════════════════════════════════════════"
echo "Launching viewer..."
echo "═══════════════════════════════════════════════════════════"
echo ""

# Launch viewer
$PYTHON_CMD "$VIEWER_PATH" "$CACHE_DIR"