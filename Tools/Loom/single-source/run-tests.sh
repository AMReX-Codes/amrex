#!/usr/bin/env bash
# Run the three single-GPU tests, then check CNS against cns-expected.txt.
set -uo pipefail
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$here/env.sh"
cd "$AMREX_BUILD"
ctest -R '^GPU_(ArenaLifetime|CNS_Sod|ExternalStream)_3d$' --output-on-failure
rc=$?
cd "$AMREX_BUILD/Tests/GPU/CNS/Sod_3d"
./Test_GPU_CNS_Sod_3d inputs-ci > cns-run.log 2>&1 || { echo "CNS run failed:"; tail -20 cns-run.log; exit 1; }
python3 "$here/check-cns.py" cns-run.log "$here/cns-expected.txt" || rc=1
exit $rc
