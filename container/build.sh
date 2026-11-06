#!/bin/bash
# Build the sRNA-IMP image from this checkout: container/build.sh [output.sif]
#   MOUNTS="/work/projects /scratch/users"   top folders of the site's file systems, created in the image
#                                            (needed where Apptainer has no overlay, e.g. ULHPC)
#   BUILD_ARGS="--fakeroot"                  apptainer build options  [--fakeroot]
#   APPTAINER_TMPDIR                         space for the build (a few GB)  [./tmp]
set -euo pipefail
cd "$(dirname "$0")/.."
APP=$(command -v apptainer || command -v singularity || true)
[ -n "$APP" ] || { echo "build.sh: apptainer / singularity not found (module load it first)" >&2; exit 1; }
VERSION=$(git describe --tags --always --dirty 2>/dev/null || echo unknown)
OUT=${1:-srna-imp_${VERSION}.sif}
echo "$VERSION" > container/VERSION
sed -e "s#__MOUNTS__#${MOUNTS:-}#g" -e "s#__VERSION__#$VERSION#g" container/srna-imp.def > container/.build.def
export APPTAINER_TMPDIR=${APPTAINER_TMPDIR:-$PWD/tmp}; mkdir -p "$APPTAINER_TMPDIR"
echo "building $OUT (version $VERSION, mount points: ${MOUNTS:-none}) ..."
$APP build ${BUILD_ARGS:---fakeroot} "$OUT" container/.build.def
rm -f container/.build.def
echo "built $(readlink -f "$OUT")"
