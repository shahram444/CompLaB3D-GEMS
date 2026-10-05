#!/usr/bin/env bash
# ---------------------------------------------------------------------------
#  Opens the .h5 viewer on macOS or Linux.
#
#      ./h5_viewer.sh                 open the window empty
#      ./h5_viewer.sh some_file.h5    open that file straight away
# ---------------------------------------------------------------------------
set -e
cd "$(dirname "$0")"

PY=${PYTHON:-python3}

if ! "$PY" -c "import h5py, numpy, matplotlib" 2>/dev/null; then
    echo "Installing what the viewer needs. This happens once."
    "$PY" -m pip install --user -r requirements.txt
fi

if ! "$PY" -c "import tkinter" 2>/dev/null; then
    echo "tkinter is missing from this Python."
    echo "  macOS, Homebrew:  brew install python-tk"
    echo "  Debian, Ubuntu:   sudo apt install python3-tk"
    exit 1
fi

exec "$PY" prt_h5_viewer.py "$@"
