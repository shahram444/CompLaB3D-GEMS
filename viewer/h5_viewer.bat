@echo off
REM ---------------------------------------------------------------------------
REM  Opens the .h5 viewer on Windows. Double click this file.
REM
REM  If a .h5 is dropped onto this file, that one opens straight away. With no
REM  file, the window opens empty and offers "Open a .h5 file".
REM ---------------------------------------------------------------------------
setlocal
cd /d "%~dp0"

where py >nul 2>nul
if %errorlevel%==0 (set PY=py -3) else (set PY=python)

%PY% -c "import h5py, numpy, matplotlib" >nul 2>nul
if not %errorlevel%==0 (
    echo Installing what the viewer needs. This happens once.
    %PY% -m pip install --user -r requirements.txt
)

%PY% prt_h5_viewer.py %1
if not %errorlevel%==0 pause
endlocal
