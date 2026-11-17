#!/bin/sh
# Compile the ModelDB model with the startup pointer fix and the shared recorder.
set -eu
variant=$1
compiler=$2
shift 2
case "$variant" in uniform|smallworld|scalefree) ;; *) exit 2 ;; esac
root=$(cd "$(dirname "$0")/.." && pwd)
base=8d5c8350fd6b602b6eea15b24088d9f67aad067e
stage=$(mktemp -d "${TMPDIR:-/tmp}/macgregor-reference.XXXXXXXX")
trap 'rm -rf "$stage"' EXIT
cd "$stage"
git -C "$root" archive -o "$stage/source.tar" "$base" \
    "$variant/main.cpp" "$variant/layer.cpp" "$variant/layer.h" "$variant/neuron.cpp" "$variant/neuron.h"
tar -xf source.tar
git apply --whitespace=nowarn --include="$variant/*" "$root/reference/recording.patch"
cp "$root"/output_options.h "$root"/run_options.h "$root"/neuron_output.h \
   "$root"/run_output.h "$root"/snapshot.h "$root"/simulation.cpp .
"$compiler" "$@" -DMACGREGOR_REFERENCE -I"$stage/$variant" \
    "$variant/main.cpp" "$variant/layer.cpp" "$variant/neuron.cpp" simulation.cpp -o macgregor-reference
mv macgregor-reference "$root/$variant/macgregor-reference"
