"""
Thin wrapper: run experiment_1.py logic on CPU with a distinct output label.
This isolates the GPU-vs-CPU effect on the same x86 machine.
"""
import subprocess, sys, os, json, shutil
from pathlib import Path

OUTPUTS = Path(__file__).parent / "experiment_1" / "outputs"

# Temporarily patch detect_platform_label to return "h_nvidia_cpu"
# by running experiment_1.py with an env var override
env = os.environ.copy()
env["PLATFORM_OVERRIDE"] = "h_nvidia_cpu"

# We need to patch the script. Simplest: import and monkey-patch, then call main.
sys.argv = ["experiment_1.py", "--device", "cpu"]

# Monkey-patch before importing
import experiment_1
experiment_1.detect_platform_label = lambda: "h_nvidia_cpu"

experiment_1.main()
