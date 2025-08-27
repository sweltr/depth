#!/usr/bin/env python3
import argparse
import glob
import difflib
import os
import shutil
import tempfile

def main():
    parser = argparse.ArgumentParser(description="Search and replace text in files")
    parser.add_argument("-i", "--input", required=True,
                        help="Input file pattern (e.g., '*.txt')")
    parser.add_argument("-f", "--find", required=True,
                        help="String to find")
    parser.add_argument("-r", "--replace", required=True,
                        help="Replacement string")
    args = parser.parse_args()

    files = glob.glob(args.input)
    if not files:
        print(f"No files found matching pattern: {args.input}")
        return

    for file in files:
        with open(file, "r", encoding="utf-8") as f:
            original = f.readlines()

        replaced = [line.replace(args.find, args.replace) for line in original]

        if original == replaced:
            continue  # skip unchanged files

        # Show only changed lines
        diff = difflib.unified_diff(
            original,
            replaced,
            fromfile=f"{file} (original)",
            tofile=f"{file} (modified)",
            lineterm="",
            n=0  # no extra context lines
        )
        print(f"\nChanges in {file}:")
        print("\n".join(diff))

        ans = input(f"Save changes to {file}? [y/N] ").strip().lower()
        if ans == "y":
            tmp_fd, tmp_path = tempfile.mkstemp()
            with os.fdopen(tmp_fd, "w", encoding="utf-8") as tmpfile:
                tmpfile.writelines(replaced)
            shutil.move(tmp_path, file)
            print(f"✔ Saved {file}")
        else:
            print(f"✘ Skipped {file}")

if __name__ == "__main__":
    main()
