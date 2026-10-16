# GeoClaw 実行ファイル(xgeoclaw)のビルド専用。実行は scripts/04_run.sh から。
CLAWMAKE = $(CLAW)/clawutil/src/Makefile.common
CLAW_PKG = geoclaw
EXE = xgeoclaw
SETRUN_FILE = setrun.py
OUTDIR = _output
SETPLOT_FILE = setplot.py
PLOTDIR = _plots
FFLAGS ?= -O2 -fopenmp
GEOLIB = $(CLAW)/geoclaw/src/2d/shallow
include $(GEOLIB)/Makefile.geoclaw
SOURCES = \
  $(CLAW)/riemann/src/rpn2_geoclaw.f \
  $(CLAW)/riemann/src/rpt2_geoclaw.f \
  $(CLAW)/riemann/src/geoclaw_riemann_utils.f
include $(CLAWMAKE)
