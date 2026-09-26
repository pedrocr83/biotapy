#!/usr/bin/env bash
# Knowledge-map freshness checks.
#
#   ./stale.sh [bundle]                  Are concepts stale vs the shared trunk?
#   ./stale.sh [bundle] --touched        Does this branch change code whose
#                                        concepts it did not update?
#
# Options:
#   --against <ref>   Trunk ref. Default: origin/main, then main, then HEAD.
#
# Exit 1 if anything is reported, so both modes work as CI gates; exit 2 when
# --touched cannot diff against the trunk ref.
set -uo pipefail

BUNDLE=".knowledge"
MODE="stale"
AGAINST=""

while [ $# -gt 0 ]; do
  case "$1" in
    --touched) MODE="touched"; shift ;;
    --against) AGAINST="$2"; shift 2 ;;
    -*) echo "Unknown option: $1" >&2; exit 2 ;;
    *) BUNDLE="$1"; shift ;;
  esac
done

if [ -z "$AGAINST" ]; then
  for ref in origin/main origin/master main master; do
    if git rev-parse --verify --quiet "$ref" >/dev/null; then AGAINST="$ref"; break; fi
  done
  AGAINST="${AGAINST:-HEAD}"
fi

[ -d "$BUNDLE" ] || { echo "No bundle at $BUNDLE" >&2; exit 1; }

python3 - "$BUNDLE" "$MODE" "$AGAINST" <<'PY'
import subprocess, sys, os, re

bundle, mode, against = sys.argv[1], sys.argv[2], sys.argv[3]

def git(*args):
    return subprocess.run(["git", *args], capture_output=True, text=True)

def frontmatter(path):
    text = open(path, encoding="utf-8").read()
    if not text.startswith("---"):
        return None
    end = text.find("\n---", 3)
    return text[3:end] if end != -1 else None

def glob_to_re(g):
    out, i = "", 0
    while i < len(g):
        c = g[i]
        if g.startswith("**/", i):   out += "(?:.*/)?"; i += 3
        elif g.startswith("**", i):  out += ".*";       i += 2
        elif c == "*":               out += "[^/]*";    i += 1
        elif c == "?":               out += "[^/]";     i += 1
        else:                        out += re.escape(c); i += 1
    return re.compile("^" + out + "$")

concepts = []
skipped = []
for root, _, files in os.walk(bundle):
    for name in sorted(files):
        if not name.endswith(".md") or name in ("index.md", "log.md"):
            continue
        path = os.path.join(root, name)
        fm = frontmatter(path)
        if fm is None:
            skipped.append((path, "no frontmatter")); continue
        commit = re.search(r"^commit:\s*['\"]?([0-9a-fA-F]+)", fm, re.M)
        paths = re.search(r"^paths:\s*\[(.*?)\]", fm, re.M | re.S)
        if not paths:
            skipped.append((path, "no paths key")); continue
        globs = [g.strip().strip("\"'") for g in paths.group(1).split(",") if g.strip()]
        concepts.append({"path": path, "commit": commit.group(1) if commit else None, "globs": globs})

if mode == "stale":
    # Has the described code moved on the trunk since the concept was written?
    # Compared against the trunk, not local HEAD, so your own in-flight work
    # does not report as staleness in someone else's map.
    stale, ok = [], 0
    for c in concepts:
        if not c["commit"]:
            skipped.append((c["path"], "no commit key")); continue
        res = git("diff", "--stat", f"{c['commit']}..{against}", "--", *c["globs"])
        if res.returncode != 0:
            reason = res.stderr.strip().splitlines()[-1] if res.stderr else "git error"
            skipped.append((c["path"], reason))
        elif res.stdout.strip():
            stale.append((c["path"], res.stdout.strip().splitlines()[-1].strip()))
        else:
            ok += 1
    if stale:
        print(f"STALE (code changed on {against} since the concept was written):")
        for p, s in stale:
            print(f"  {p}\n      {s}")
    if skipped:
        print("\nUNCHECKABLE:")
        for p, why in skipped:
            print(f"  {p}  ({why})")
    print(f"\n{ok} current, {len(stale)} stale, {len(skipped)} uncheckable")
    sys.exit(1 if stale else 0)

# mode == "touched": which concepts cover code this branch changed but did not update?
def git_or_exit(*args):
    res = git(*args)
    if res.returncode != 0:
        print(f"cannot compare against {against}: {res.stderr.strip()}", file=sys.stderr)
        sys.exit(2)
    return res.stdout

base = git_or_exit("merge-base", against, "HEAD").strip()
changed = [f for f in git_or_exit("diff", "--name-only", f"{base}..HEAD").splitlines() if f]
edited_concepts = {f for f in changed if f.startswith(bundle.rstrip("/") + "/")}
code_changed = [f for f in changed if f not in edited_concepts]

flagged = []
for c in concepts:
    if c["path"] in edited_concepts:
        continue
    pats = [glob_to_re(g) for g in c["globs"]]
    hits = [f for f in code_changed if any(p.match(f) for p in pats)]
    if hits:
        flagged.append((c["path"], hits))

if flagged:
    print(f"This branch changes code covered by concepts it did not update (base {base[:9]}):\n")
    for p, hits in flagged:
        print(f"  {p}")
        for h in hits[:5]:
            print(f"      {h}")
        if len(hits) > 5:
            print(f"      ... and {len(hits) - 5} more")
    print("\nUpdate the concept if the change invalidates something it states,")
    print("and refresh its `commit`. Most changes will not: implementation detail")
    print("does not belong in the map. Confirm rather than reflexively editing.")
sys.exit(1 if flagged else 0)
PY
