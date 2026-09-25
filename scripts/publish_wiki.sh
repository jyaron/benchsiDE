#!/usr/bin/env bash
# Publish wiki/ to the GitHub repository wiki (a separate git repository, <repo>.wiki.git).
# Prerequisite: the wiki must exist. Create it once in the browser (repository > Wiki >
# "Create the first page" > Save); that creates <repo>.wiki.git.
# Usage: scripts/publish_wiki.sh [owner/repo]     (default: jyaron/benchsiDE)
set -euo pipefail
REPO="${1:-jyaron/benchsiDE}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

git clone --quiet "https://github.com/${REPO}.wiki.git" "$WORK/wiki"
# replace wiki contents with the versioned source, rewriting repository-relative links:
#   (Page.md)            -> (Page)                          wiki page links
#   (../docs/x.md)       -> (https://github.com/REPO/blob/main/docs/x.md)
#   (../VALIDATION.md)   -> (https://github.com/REPO/blob/main/VALIDATION.md)
find "$WORK/wiki" -maxdepth 1 -name '*.md' -delete
for f in "$ROOT"/wiki/*.md; do
  sed -E \
    -e "s#\]\(\.\./([^)]+)\)#](https://github.com/${REPO}/blob/main/\1)#g" \
    -e 's#\]\(([A-Za-z0-9_-]+)\.md(\#[^)]*)?\)#](\1\2)#g' \
    "$f" > "$WORK/wiki/$(basename "$f")"
done
cd "$WORK/wiki"
git add -A
if git diff --cached --quiet; then
  echo "wiki already up to date"
else
  git commit --quiet -m "Sync knowledgebase from main ($(cd "$ROOT" && git rev-parse --short HEAD 2>/dev/null || echo local))"
  git push --quiet
  echo "wiki published: https://github.com/${REPO}/wiki"
fi
