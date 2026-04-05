"""Generate hardware manifest for the current platform (auto-detected)."""
import json
import platform
import subprocess
import datetime
from pathlib import Path


def detect_platform_label():
    """Auto-detect platform label based on available hardware."""
    import torch
    if torch.cuda.is_available():
        return "H_nvidia"
    elif hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        return "H_mac"
    else:
        return "H_cpu"


def generate_manifest():
    import torch
    import transformers
    import numpy as np

    label = detect_platform_label()

    manifest = {
        "platform_label": label,
        "timestamp": datetime.datetime.now().isoformat(),
        "os": platform.platform(),
        "os_version": platform.version(),
        "processor": platform.processor(),
        "python_version": platform.python_version(),
        "pytorch_version": torch.__version__,
        "transformers_version": transformers.__version__,
        "numpy_version": np.__version__,
        "cuda_available": torch.cuda.is_available(),
        "mps_available": hasattr(torch.backends, "mps") and torch.backends.mps.is_available(),
    }

    if torch.cuda.is_available():
        props = torch.cuda.get_device_properties(0)
        manifest.update({
            "gpu_name": props.name,
            "gpu_memory_mb": props.total_memory // (1024 * 1024),
            "gpu_compute_capability": f"{props.major}.{props.minor}",
            "cuda_version": torch.version.cuda,
            "cudnn_version": str(torch.backends.cudnn.version()),
            "bf16_supported": torch.cuda.is_bf16_supported(),
        })
        try:
            result = subprocess.run(
                ["nvidia-smi", "--query-gpu=driver_version", "--format=csv,noheader"],
                capture_output=True, text=True, timeout=5,
            )
            if result.returncode == 0:
                manifest["nvidia_driver_version"] = result.stdout.strip()
        except Exception:
            pass

    elif manifest["mps_available"]:
        # Apple Silicon info
        try:
            result = subprocess.run(
                ["sysctl", "-n", "machdep.cpu.brand_string"],
                capture_output=True, text=True, timeout=5,
            )
            if result.returncode == 0:
                manifest["chip"] = result.stdout.strip()
        except Exception:
            pass
        try:
            result = subprocess.run(
                ["system_profiler", "SPHardwareDataType"],
                capture_output=True, text=True, timeout=10,
            )
            if result.returncode == 0:
                for line in result.stdout.splitlines():
                    line = line.strip()
                    if line.startswith("Memory:"):
                        manifest["memory"] = line.split(":", 1)[1].strip()
        except Exception:
            pass
        manifest["bf16_supported"] = True  # All Apple Silicon supports BF16

    return manifest


if __name__ == "__main__":
    m = generate_manifest()
    label = m["platform_label"].lower()
    out = Path(__file__).parent / f"hardware_manifest_{label}.json"
    out.write_text(json.dumps(m, indent=2) + "\n")
    print(json.dumps(m, indent=2))
    print(f"\nSaved to {out}")
