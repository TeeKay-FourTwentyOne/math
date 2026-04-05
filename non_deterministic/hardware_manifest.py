"""Generate hardware manifest for the H_nvidia platform."""
import json
import platform
import subprocess
import datetime
from pathlib import Path


def generate_manifest():
    import torch
    import transformers
    import numpy as np

    manifest = {
        "platform_label": "H_nvidia",
        "timestamp": datetime.datetime.now().isoformat(),
        "os": platform.platform(),
        "os_version": platform.version(),
        "processor": platform.processor(),
        "python_version": platform.python_version(),
        "pytorch_version": torch.__version__,
        "transformers_version": transformers.__version__,
        "numpy_version": np.__version__,
        "cuda_available": torch.cuda.is_available(),
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

    return manifest


if __name__ == "__main__":
    m = generate_manifest()
    out = Path(__file__).parent / "hardware_manifest.json"
    out.write_text(json.dumps(m, indent=2))
    print(json.dumps(m, indent=2))
    print(f"\nSaved to {out}")
