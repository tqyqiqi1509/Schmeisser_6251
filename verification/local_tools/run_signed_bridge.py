#!/usr/bin/env python3
"""Generate the signed scalar certificate, then verify the actual W8 bridge.

The scalar receipt is generated in this invocation. No stored acceptance
receipt is a mathematical input.
"""
from pathlib import Path
import subprocess
import sys

verification = Path(__file__).resolve().parent.parent
source = verification / "source"
for script in (
    "computations/degree6251/certify_signed_theta_q06.py",
    "agents/small12501_general/review_tang_W8_signed_field_independent.py",
):
    subprocess.run([sys.executable, str(source / script)], cwd=source, check=True)
