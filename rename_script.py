"""
Rename epimatch -> epimatch across the repo:
- epimatch/ package directory -> epimatch/
- "epimatch" (lowercase, module/package identifier) -> "epimatch"
- "EpiMatch" (display/class name, e.g. EpiMatchDataset) -> "EpiMatch"
- GitHub repo URL (gsyadav2000/EpiMatch) is left untouched - that's a
  separate decision, see the instructions after running this.
"""

import os
import shutil

ROOT = "."
GITHUB_URL_MARKER = "github.com/gsyadav2000/AtlasX"
PLACEHOLDER = "github.com/gsyadav2000/AtlasX"

if os.path.isdir(os.path.join(ROOT, "epimatch")) and not os.path.isdir(os.path.join(ROOT, "epimatch")):
    shutil.move(os.path.join(ROOT, "epimatch"), os.path.join(ROOT, "epimatch"))
    print("Renamed directory: epimatch/ -> epimatch/")
else:
    print("Skipped directory rename (already done or epimatch/ not found)")

extensions = (".py", ".toml", ".md", ".txt")
changed_files = []

for dirpath, dirnames, filenames in os.walk(ROOT):
    dirnames[:] = [d for d in dirnames if d not in (".git", "__pycache__", ".venv", "venv")]

    for filename in filenames:
        if not filename.endswith(extensions):
            continue

        filepath = os.path.join(dirpath, filename)

        with open(filepath, "r", encoding="utf-8") as f:
            original = f.read()

        text = original
        text = text.replace(GITHUB_URL_MARKER, PLACEHOLDER)
        text = text.replace("epimatch", "epimatch")
        text = text.replace("EpiMatch", "EpiMatch")
        text = text.replace(PLACEHOLDER, GITHUB_URL_MARKER)

        if text != original:
            with open(filepath, "w", encoding="utf-8") as f:
                f.write(text)
            changed_files.append(filepath)

print(f"\nModified {len(changed_files)} files.")