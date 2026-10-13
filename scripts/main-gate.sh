#!/bin/sh
# Refuse a push that moves main unless the commit is a checked release.
# Runs as a pre-push hook; CLAUDE.md says how to register it. Every other
# push passes through untouched.

set -u

fail() {
    echo "main-gate: $*" >&2
    exit 1
}

# git hands a pre-push hook one line per ref:
#   <local ref> <local sha> <remote ref> <remote sha>
main_sha=
while read -r _ local_sha remote_ref _; do
    if [ "$remote_ref" = refs/heads/main ]; then
        main_sha=$local_sha
    fi
done

[ -n "$main_sha" ] || exit 0

case $main_sha in
*[!0]*) ;;
*) fail "refusing to delete main" ;;
esac

for tool in git cargo cargo-semver-checks; do
    command -v "$tool" >/dev/null 2>&1 ||
        fail "$tool is not installed, and the checks on main need it"
done

# the checks run on the pushed commit itself, in a worktree of its own,
# so neither the branch checked out nor uncommitted edits can sway them.
# its build directory stays between pushes, so a check is incremental
top=$(git rev-parse --show-toplevel)

# git may export these to a hook, and they would point the commands below
# back at the checkout the push came from
unset GIT_DIR GIT_WORK_TREE GIT_INDEX_FILE

gate="$top/target/main-gate"
tree="$gate/tree"

git worktree remove --force "$tree" 2>/dev/null
git worktree add --quiet --detach "$tree" "$main_sha" ||
    fail "could not check out $main_sha to check it"
trap 'git -C "$top" worktree remove --force "$tree"' EXIT
cd "$tree" || fail "could not enter $tree"
export CARGO_TARGET_DIR="$gate/target"

tags=$(git tag --points-at "$main_sha" | grep -E '^(libsail|sail)-v[0-9]')
[ -n "$tags" ] ||
    fail "main moves only to a release; tag the commit <crate>-v<version>"

for tag in $tags; do
    crate=${tag%%-v*}
    version=${tag#*-v}

    [ "$(git cat-file -t "$tag")" = tag ] ||
        fail "$tag is a lightweight tag; a release tag is annotated"
    grep -qx "version = \"$version\"" "$crate/Cargo.toml" ||
        fail "$tag does not match the version in $crate/Cargo.toml"
    grep -qF "## [$version]" "$crate/CHANGELOG.md" ||
        fail "$crate/CHANGELOG.md has no heading for $version"

    if [ "$crate" = libsail ]; then
        previous=$(git tag --list 'libsail-v*' --sort=-v:refname |
            grep -vx "$tag" | head -n 1)
        [ -n "$previous" ] ||
            fail "no earlier libsail tag to compare $tag against"
        cargo semver-checks -p libsail --baseline-rev "$previous" ||
            fail "cargo semver-checks rejects $tag against $previous"
    fi
done

cargo fmt --all --check ||
    fail "cargo fmt --all --check failed"
cargo clippy --workspace --all-targets --all-features -- -D warnings ||
    fail "cargo clippy reported warnings or errors"
cargo test --workspace --all-features ||
    fail "cargo test failed"
